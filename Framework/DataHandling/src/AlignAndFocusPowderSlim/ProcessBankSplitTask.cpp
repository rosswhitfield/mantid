// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2025 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "MantidDataHandling/AlignAndFocusPowderSlim/ProcessBankSplitTask.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/ProcessEventsTask.h"
#include "MantidKernel/Logger.h"
#include "MantidKernel/ParallelMinMax.h"
#include "MantidKernel/Timer.h"
#include "MantidKernel/Unit.h"
#include "MantidNexus/H5Util.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

#include <ranges>

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

namespace {

const std::string MICROSEC("microseconds");

// Logger for this class
auto g_log = Kernel::Logger("ProcessBankSplitTask");

} // namespace
ProcessBankSplitTask::ProcessBankSplitTask(
    std::vector<std::string> &bankEntryNames, H5::H5File &h5file, const bool is_time_filtered,
    std::vector<int> &workspaceIndices, std::vector<API::MatrixWorkspace_sptr> &wksps,
    const std::map<detid_t, double> &calibration, const std::map<detid_t, double> &scale_at_sample,
    const std::set<detid_t> &masked, const size_t events_per_chunk, const size_t grainsize_event,
    std::vector<std::pair<int, PulseROI>> target_to_pulse_indices, std::shared_ptr<API::Progress> &progress)
    : ProcessBankTaskBase(h5file, bankEntryNames, calibration, scale_at_sample, masked, events_per_chunk,
                          grainsize_event, progress),
      m_loader(is_time_filtered, {}, target_to_pulse_indices), m_workspaceIndices(workspaceIndices), m_wksps(wksps) {}

void ProcessBankSplitTask::operator()(const tbb::blocked_range<size_t> &range) const {
  auto entry = m_h5file.openGroup("entry"); // type=NXentry
  for (size_t wksp_index = range.begin(); wksp_index < range.end(); ++wksp_index) {
    const auto &bankName = m_bankEntries[wksp_index];
    // empty bank names indicate spectra to skip; control should never get here, but just in case
    if (bankName.empty()) {
      continue;
    }
    Kernel::Timer timer;
    g_log.debug() << bankName << " start" << std::endl;

    // open the bank
    auto event_group = entry.openGroup(bankName); // type=NXevent_data

    // skip empty dataset
    auto tof_SDS = event_group.openDataSet(NxsFieldNames::TIME_OF_FLIGHT);
    const int64_t total_events = static_cast<size_t>(tof_SDS.getSpace().getSelectNpoints());
    if (total_events == 0) {
      m_progress->report();
      continue;
    }

    auto eventSplitRanges = m_loader.getEventIndexSplitRanges(event_group, total_events);

    // Get all spectra for this bank.
    // Create temporary y arrays for each workspace.
    std::vector<API::ISpectrum *> spectra;
    std::vector<std::vector<uint32_t>> y_temps;
    for (const auto &wksp : m_wksps) {
      spectra.push_back(&wksp->getSpectrum(wksp_index));
      y_temps.emplace_back(spectra.back()->dataY().size());
    }

    // create object so bank calibration can be re-used
    std::unique_ptr<BankCalibration> calibration = nullptr;

    // get handle to the data
    auto detID_SDS = event_group.openDataSet(NxsFieldNames::DETID);
    // auto tof_SDS = event_group.openDataSet(NxsFieldNames::TIME_OF_FLIGHT);
    // and the units
    std::string tof_unit;
    Nexus::H5Util::readStringAttribute(tof_SDS, "units", tof_unit);
    const double time_conversion = Kernel::Units::timeConversionValue(tof_unit, MICROSEC);

    // declare arrays once so memory can be reused
    auto event_detid = std::make_unique<std::vector<uint32_t>>();       // uint32 for ORNL nexus file
    auto event_time_of_flight = std::make_unique<std::vector<float>>(); // float for ORNL nexus files

    // read parts of the bank at a time until all events are processed
    while (!eventSplitRanges.empty()) {
      // Create offsets and slab sizes for the next chunk of events.
      auto chunk = prepareNextSplitChunk(eventSplitRanges);

      // log the event ranges being processed
      logChunkInfo(bankName, chunk);

      // load detid and tof at the same time
      tbb::parallel_invoke(
          [&] { // load detid
            // event_detid->clear();
            m_loader.loadData(detID_SDS, event_detid, chunk.offsets, chunk.slabsizes);
            // immediately find min/max to allow for other things to read disk
            const auto [minval, maxval] = Mantid::Kernel::parallel_minmax(event_detid, m_grainsize_event);
            // only recreate calibration if it doesn't already have the useful information
            updateCalibration(calibration, minval, maxval, time_conversion);
          },
          [&] { // load time-of-flight
            // event_time_of_flight->clear();
            m_loader.loadData(tof_SDS, event_time_of_flight, chunk.offsets, chunk.slabsizes);
          });

      // loop over targets
      tbb::parallel_for(
          tbb::blocked_range<size_t>(0, m_workspaceIndices.size()), [&](const tbb::blocked_range<size_t> &r) {
            for (size_t idx = r.begin(); idx != r.end(); ++idx) {
              int i = m_workspaceIndices[idx];

              // Precompute indices for this target
              std::vector<size_t> indices;
              for (const auto &pair : chunk.relative_target_ranges) {
                if (pair.first == static_cast<int>(i)) {
                  auto [start, end] = pair.second;
                  for (size_t k = start; k < end; ++k) {
                    indices.push_back(k);
                  }
                }
              }

              auto event_id_view_for_target =
                  indices | std::views::transform([&event_detid](const auto &k) { return (*event_detid)[k]; });
              auto event_tof_view_for_target = indices | std::views::transform([&event_time_of_flight](const auto &k) {
                                                 return (*event_time_of_flight)[k];
                                               });

              ProcessEventsTask task(&event_id_view_for_target, &event_tof_view_for_target, calibration.get(),
                                     &spectra[idx]->readX());

              const tbb::blocked_range<size_t> range_info(0, indices.size(), m_grainsize_event);
              tbb::parallel_reduce(range_info, task);

              std::transform(y_temps[idx].begin(), y_temps[idx].end(), task.y_temp.begin(), y_temps[idx].begin(),
                             std::plus<uint32_t>());
            }
          });
    }

    // copy the data out into the correct spectrum and calculate errors
    tbb::parallel_for(size_t(0), m_wksps.size(),
                      [&](size_t i) { copyResultsAndCalculateErrors(*spectra[i], y_temps[i]); });

    g_log.debug() << bankName << " stop " << timer << std::endl;
    m_progress->report();
  }
}
}; // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
