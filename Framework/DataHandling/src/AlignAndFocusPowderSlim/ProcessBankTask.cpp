// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2025 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "ProcessBankTask.h"
#include "MantidKernel/Logger.h"
#include "MantidKernel/Timer.h"
#include "MantidKernel/Unit.h"
#include "MantidNexus/H5Util.h"
#include "ProcessEventsTask.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

namespace {

const std::string MICROSEC("microseconds");

// Logger for this class
auto g_log = Kernel::Logger("ProcessBankTask");

template <typename Type> class MinMax {
  const std::vector<Type> *vec;

public:
  Type minval;
  Type maxval;
  void operator()(const tbb::blocked_range<size_t> &range) {
    const auto [minele, maxele] = std::minmax_element(vec->cbegin() + range.begin(), vec->cbegin() + range.end());
    if (*minele < minval)
      minval = *minele;
    if (*maxele > maxval)
      maxval = *maxele;
  }

  // copy min/max from the other. we're all friends
  MinMax(MinMax &other, tbb::split) : vec(other.vec), minval(other.minval), maxval(other.maxval) {}

  // set the min=max=first element supplied
  MinMax(const std::vector<Type> *vec) : vec(vec), minval(vec->front()), maxval(vec->front()) {}

  void join(const MinMax &other) {
    if (other.minval < minval)
      minval = other.minval;
    if (other.maxval > maxval)
      maxval = other.maxval;
  }
};

template <typename Type> std::pair<Type, Type> parallel_minmax(const std::vector<Type> *vec, const size_t grainsize) {
  if (vec->size() < grainsize) {
    const auto [minval, maxval] = std::minmax_element(vec->cbegin(), vec->cend());
    return std::make_pair(*minval, *maxval);
  } else {
    MinMax<Type> finder(vec);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, vec->size(), grainsize), finder);
    return std::make_pair(finder.minval, finder.maxval);
  }
}

} // namespace
ProcessBankTask::ProcessBankTask(std::vector<std::string> &bankEntryNames, H5::H5File &h5file,
                                 const bool is_time_filtered, API::MatrixWorkspace_sptr &wksp,
                                 const std::map<detid_t, double> &calibration, const std::set<detid_t> &masked,
                                 const size_t events_per_chunk, const size_t grainsize_event,
                                 std::vector<PulseROI> pulse_indices, std::shared_ptr<API::Progress> &progress)
    : m_h5file(h5file), m_bankEntries(bankEntryNames), m_loader(is_time_filtered, pulse_indices), m_wksp(wksp),
      m_calibration(calibration), m_masked(masked), m_events_per_chunk(events_per_chunk),
      m_grainsize_event(grainsize_event), m_progress(progress) {}

void ProcessBankTask::operator()(const tbb::blocked_range<size_t> &range) const {
  auto entry = m_h5file.openGroup("entry"); // type=NXentry
  for (size_t wksp_index = range.begin(); wksp_index < range.end(); ++wksp_index) {
    const auto &bankName = m_bankEntries[wksp_index];
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
    std::vector<std::atomic_uint32_t> y_temp(spectrum.dataY().size());
    // std::vector<uint32_t> y_temp(spectrum.dataY().size());

    // task group allows for separate of disk read from processing
    // tbb::task_group_context tgroupcontext; // needed by parallel_reduce
    // tbb::task_group tgroup(tgroupcontext);

    // create object so bank calibration can be re-used

    // get handle to the data
    auto detID_SDS = event_group.openDataSet(NxsFieldNames::DETID);
    // auto tof_SDS = event_group.openDataSet(NxsFieldNames::TIME_OF_FLIGHT);
    // and the units
    std::string tof_unit;
    Nexus::H5Util::readStringAttribute(tof_SDS, "units", tof_unit);
    const double time_conversion = Kernel::Units::timeConversionValue(tof_unit, MICROSEC);

    // Prepare all chunks for this bank
    std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> chunk_offsets_slabs;
    while (!eventRanges.empty()) {
      std::vector<size_t> offsets, slabsizes;
      size_t total_events_to_read = 0;
      while (!eventRanges.empty() && total_events_to_read < m_events_per_chunk) {
        auto eventRange = eventRanges.top();
        eventRanges.pop();
        size_t range_size = eventRange.second - eventRange.first;
        size_t remaining_chunk = m_events_per_chunk - total_events_to_read;
        if (range_size > remaining_chunk) {
          offsets.push_back(eventRange.first);
          slabsizes.push_back(remaining_chunk);
          total_events_to_read += remaining_chunk;
          eventRanges.emplace(eventRange.first + remaining_chunk, eventRange.second);
          break;
        } else {
          offsets.push_back(eventRange.first);
          slabsizes.push_back(range_size);
          total_events_to_read += range_size;
        }
      }
      if (!offsets.empty())
        chunk_offsets_slabs.emplace_back(std::move(offsets), std::move(slabsizes));
    }

    // Parallelize over all chunks for this bank
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, chunk_offsets_slabs.size()), [&](const tbb::blocked_range<size_t> &chunk_range) {
          for (size_t chunk_idx = chunk_range.begin(); chunk_idx != chunk_range.end(); ++chunk_idx) {
            const auto &[offsets, slabsizes] = chunk_offsets_slabs[chunk_idx];
            std::unique_ptr<BankCalibration> calibration = nullptr;

            // Local vectors for this chunk
            auto event_detid = std::make_unique<std::vector<uint32_t>>();
            auto event_time_of_flight = std::make_unique<std::vector<float>>();

            // Load data in parallel
            tbb::parallel_invoke(
                [&] {
                  m_loader.loadData(detID_SDS, event_detid, offsets, slabsizes);
                  const auto [minval, maxval] = parallel_minmax(event_detid.get(), m_grainsize_event);
                  calibration =
                      std::make_unique<BankCalibration>(static_cast<detid_t>(minval), static_cast<detid_t>(maxval),
                                                        time_conversion, m_calibration, m_masked);
                },
                [&] { m_loader.loadData(tof_SDS, event_time_of_flight, offsets, slabsizes); });

            // Process events
            ProcessEventsTask task(event_detid.get(), event_time_of_flight.get(), calibration.get(), &spectrum.readX());
            const tbb::blocked_range<size_t> range_info(0, event_time_of_flight->size(), m_grainsize_event);
            tbb::parallel_reduce(range_info, task);

            // Atomically accumulate results into shared y_temp
            for (size_t i = 0; i < y_temp.size(); ++i) {
              y_temp[i].fetch_add(task.y_temp[i], std::memory_order_relaxed);
            }
          }
        });

    // copy the data out into the correct spectrum
    auto &y_values = spectrum.dataY();
    std::copy(y_temp.cbegin(), y_temp.cend(), y_values.begin());

    g_log.debug() << bankName << " stop " << timer << std::endl;
    m_progress->report();
  }
}
}; // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
