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
                                 const bool is_time_filtered, SpectraProcessingData &processingData,
                                 const std::map<detid_t, double> &calibration,
                                 const std::map<detid_t, double> &scale_at_sample,
                                 const std::map<size_t, std::vector<detid_t>> &grouping,
                                 const std::set<detid_t> &masked, const size_t events_per_chunk,
                                 const size_t grainsize_event, std::vector<PulseROI> pulse_indices,
                                 std::shared_ptr<API::Progress> &progress)
    : m_h5file(h5file), m_bankEntries(bankEntryNames), m_loader(is_time_filtered, pulse_indices),
      m_processingData(processingData), m_calibration(calibration), m_scale_at_sample(scale_at_sample),
      m_grouping(grouping), m_masked(masked), m_events_per_chunk(events_per_chunk), m_grainsize_event(grainsize_event),
      m_progress(progress) {}

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
      // This will read at most m_events_per_chunk events from the file
      // and will split the ranges if necessary for the next iteration.
      std::vector<size_t> offsets;
      std::vector<size_t> slabsizes;

      size_t total_events_to_read = 0;
      // Process the event ranges until we reach the desired number of events to read or run out of ranges
      while (!eventRanges.empty() && total_events_to_read < m_events_per_chunk) {
        // Get the next event range from the stack
        auto eventRange = eventRanges.top();
        eventRanges.pop();

        size_t range_size = eventRange.second - eventRange.first;
        size_t remaining_chunk = m_events_per_chunk - total_events_to_read;

        // If the range size is larger than the remaining chunk, we need to split it
        if (range_size > remaining_chunk) {
          // Split the range: process only part of it now, push the rest back for later
          offsets.push_back(eventRange.first);
          slabsizes.push_back(remaining_chunk);
          total_events_to_read += remaining_chunk;
          // Push the remainder of the range back to the front for next iteration
          eventRanges.emplace(eventRange.first + remaining_chunk, eventRange.second);
          break;
        } else {
          offsets.push_back(eventRange.first);
          slabsizes.push_back(range_size);
          total_events_to_read += range_size;
          // Continue to next range
        }
      }

      // log the event ranges being processed
      std::ostringstream oss;
      oss << "Processing " << bankName << " with " << total_events_to_read << " events in the ranges: ";
      for (size_t i = 0; i < offsets.size(); ++i) {
        oss << "[" << offsets[i] << ", " << (offsets[i] + slabsizes[i]) << "), ";
      }
      oss << "\n";
      g_log.debug() << oss.str();

      detid_t minval, maxval;

      // load detid and tof at the same time
      tbb::parallel_invoke(
          [&] { // load detid
            // event_detid->clear();
            m_loader.loadData(detID_SDS, event_detid, offsets, slabsizes);
            // immediately find min/max to allow for other things to read disk
            std::tie(minval, maxval) = Mantid::Kernel::parallel_minmax(event_detid, m_grainsize_event);

          },
          [&] { // load time-of-flight
            // event_time_of_flight->clear();
            m_loader.loadData(tof_SDS, event_time_of_flight, offsets, slabsizes);
          });

      // Loop over all groups in m_grouping for this bank
      for (const auto &[group_index, detids_in_group] : m_grouping) {
        calibration = std::make_unique<BankCalibration>(static_cast<detid_t>(minval), static_cast<detid_t>(maxval),
                                                        time_conversion, detids_in_group, m_calibration,
                                                        m_scale_at_sample, m_masked);

        // Create a local task for this thread
        ProcessEventsTask task(event_detid.get(), event_time_of_flight.get(), calibration.get(),
                               m_processingData.binedges[group_index]);

        // Non-blocking processing of the events
        const tbb::blocked_range<size_t> range_info(0, event_time_of_flight->size(), m_grainsize_event);
        tbb::parallel_reduce(range_info, task);

        // Accumulate results into shared y_temp to combine local histograms
        // Use atomic fetch_add to accumulate results into shared y_temp
        for (size_t i = 0; i < m_processingData.counts[group_index].size(); ++i) {
          m_processingData.counts[group_index][i].fetch_add(task.y_temp[i], std::memory_order_relaxed);
        }
      }
    }

    g_log.debug() << bankName << " stop " << timer << std::endl;
    m_progress->report();
  }
}
}; // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
