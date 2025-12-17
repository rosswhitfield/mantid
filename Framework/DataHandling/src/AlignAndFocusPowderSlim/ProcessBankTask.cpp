// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2025 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "MantidDataHandling/AlignAndFocusPowderSlim/ProcessBankTask.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/ProcessEventsTask.h"
#include "MantidKernel/Logger.h"
#include "MantidKernel/ParallelMinMax.h"
#include "MantidKernel/Timer.h"
#include "MantidKernel/Unit.h"
#include "MantidNexus/H5Util.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

namespace {

const std::string MICROSEC("microseconds");

// Logger for this class
auto g_log = Kernel::Logger("ProcessBankTask");

} // namespace
ProcessBankTask::ProcessBankTask(std::vector<std::string> &bankEntryNames, H5::H5File &h5file,
                                 const bool is_time_filtered, API::MatrixWorkspace_sptr &wksp,
                                 const std::map<detid_t, double> &calibration,
                                 const std::map<detid_t, double> &scale_at_sample, const std::set<detid_t> &masked,
                                 const size_t events_per_chunk, const size_t grainsize_event,
                                 std::vector<PulseROI> pulse_indices, std::shared_ptr<API::Progress> &progress)
    : ProcessBankTaskBase(h5file, bankEntryNames, calibration, scale_at_sample, masked, events_per_chunk,
                          grainsize_event, progress),
      m_loader(is_time_filtered, pulse_indices), m_wksp(wksp) {}

void ProcessBankTask::operator()(const tbb::blocked_range<size_t> &range) const {
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

    auto eventRanges = m_loader.getEventIndexRanges(event_group, total_events);

    // create a histogrammer to process the events
    auto &spectrum = m_wksp->getSpectrum(wksp_index);

    // std::atomic allows for multi-threaded accumulation and who cares about floats when you are just
    // counting things
    // std::vector<std::atomic_uint32_t> y_temp(spectrum.dataY().size())
    std::vector<uint32_t> y_temp(spectrum.dataY().size());

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
    while (!eventRanges.empty()) {
      // Create offsets and slab sizes for the next chunk of events.
      auto chunk = prepareNextChunk(eventRanges);

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

      // Create a local task for this thread
      ProcessEventsTask task(event_detid.get(), event_time_of_flight.get(), calibration.get(), &spectrum.readX());

      // Non-blocking processing of the events
      const tbb::blocked_range<size_t> range_info(0, event_time_of_flight->size(), m_grainsize_event);
      tbb::parallel_reduce(range_info, task);

      // Accumulate results into shared y_temp to combine local histograms
      std::transform(y_temp.begin(), y_temp.end(), task.y_temp.begin(), y_temp.begin(), std::plus<uint32_t>());
    }

    // copy the data out into the correct spectrum and calculate errors
    copyResultsAndCalculateErrors(spectrum, y_temp);

    g_log.debug() << bankName << " stop " << timer << std::endl;
    m_progress->report();
  }
}
}; // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
