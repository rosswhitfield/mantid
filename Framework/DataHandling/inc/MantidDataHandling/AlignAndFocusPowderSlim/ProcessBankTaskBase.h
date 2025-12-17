// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2025 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/MatrixWorkspace.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/BankCalibration.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/NexusLoader.h"
#include "MantidGeometry/IDTypes.h"
#include "MantidKernel/Logger.h"
#include <H5Cpp.h>
#include <MantidAPI/Progress.h>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stack>
#include <vector>

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

/// Base class containing common bank processing functionality
class ProcessBankTaskBase {
protected:
  ProcessBankTaskBase(H5::H5File &h5file, const std::vector<std::string> &bankEntries,
                      const std::map<detid_t, double> &calibration, const std::map<detid_t, double> &scale_at_sample,
                      const std::set<detid_t> &masked, const size_t events_per_chunk, const size_t grainsize_event,
                      std::shared_ptr<API::Progress> &progress)
      : m_h5file(h5file), m_bankEntries(bankEntries), m_calibration(calibration), m_scale_at_sample(scale_at_sample),
        m_masked(masked), m_events_per_chunk(events_per_chunk), m_grainsize_event(grainsize_event),
        m_progress(progress) {}

  /// Common structure for chunking information
  struct ChunkInfo {
    std::vector<size_t> offsets;
    std::vector<size_t> slabsizes;
    size_t total_events_to_read;
  };

  /// Structure for split chunking with target information
  struct SplitChunkInfo : ChunkInfo {
    std::vector<std::pair<int, EventROI>> relative_target_ranges;
  };

  /// Prepare next chunk of events to read from disk
  template <typename RangeStack> ChunkInfo prepareNextChunk(RangeStack &ranges) const {
    ChunkInfo chunk;
    chunk.total_events_to_read = 0;

    // Process the event ranges until we reach the desired number of events to read or run out of ranges
    while (!ranges.empty() && chunk.total_events_to_read < m_events_per_chunk) {
      auto eventRange = ranges.top();
      ranges.pop();

      size_t range_size = eventRange.second - eventRange.first;
      size_t remaining_chunk = m_events_per_chunk - chunk.total_events_to_read;

      // If the range size is larger than the remaining chunk, we need to split it
      if (range_size > remaining_chunk) {
        // Split the range: process only part of it now, push the rest back for later
        chunk.offsets.push_back(eventRange.first);
        chunk.slabsizes.push_back(remaining_chunk);
        chunk.total_events_to_read += remaining_chunk;
        // Push the remainder of the range back to the front for next iteration
        ranges.emplace(eventRange.first + remaining_chunk, eventRange.second);
        break;
      } else {
        chunk.offsets.push_back(eventRange.first);
        chunk.slabsizes.push_back(range_size);
        chunk.total_events_to_read += range_size;
        // Continue to next range
      }
    }

    return chunk;
  }

  /// Prepare next chunk with split target information
  template <typename SplitRangeStack> SplitChunkInfo prepareNextSplitChunk(SplitRangeStack &ranges) const {
    SplitChunkInfo chunk;
    chunk.total_events_to_read = 0;

    // Process the event ranges until we reach the desired number of events to read or run out of ranges
    while (!ranges.empty() && chunk.total_events_to_read < m_events_per_chunk) {
      auto [target, eventRange] = ranges.top();
      ranges.pop();

      size_t range_size = eventRange.second - eventRange.first;
      size_t remaining_chunk = m_events_per_chunk - chunk.total_events_to_read;

      // If the range size is larger than the remaining chunk, we need to split it
      if (range_size > remaining_chunk) {
        // Split the range: process only part of it now, push the rest back for later
        chunk.relative_target_ranges.emplace_back(
            target, EventROI(chunk.total_events_to_read, chunk.total_events_to_read + remaining_chunk));
        chunk.offsets.push_back(eventRange.first);
        chunk.slabsizes.push_back(remaining_chunk);
        chunk.total_events_to_read += remaining_chunk;
        // Push the remainder of the range back to the front for next iteration
        ranges.emplace(target, EventROI(eventRange.first + remaining_chunk, eventRange.second));
        break;
      } else {
        chunk.relative_target_ranges.emplace_back(
            target, EventROI(chunk.total_events_to_read, chunk.total_events_to_read + range_size));
        chunk.offsets.push_back(eventRange.first);
        chunk.slabsizes.push_back(range_size);
        chunk.total_events_to_read += range_size;
        // Continue to next range
      }
    }

    return chunk;
  }

  /// Log debug information about the chunk being processed
  void logChunkInfo(const std::string &bankName, const ChunkInfo &chunk) const {
    std::ostringstream oss;
    oss << "Processing " << bankName << " with " << chunk.total_events_to_read << " events in the ranges: ";
    for (size_t i = 0; i < chunk.offsets.size(); ++i) {
      oss << "[" << chunk.offsets[i] << ", " << (chunk.offsets[i] + chunk.slabsizes[i]) << "), ";
    }
    oss << "\n";
    g_log.debug() << oss.str();
  }

  /// Update or create calibration object if needed based on min/max detector IDs
  void updateCalibration(std::unique_ptr<BankCalibration> &calibration, uint32_t minval, uint32_t maxval,
                         double time_conversion) const {
    if ((!calibration) || (calibration->idmin() > static_cast<detid_t>(minval)) ||
        (calibration->idmax() < static_cast<detid_t>(maxval))) {
      calibration = std::make_unique<BankCalibration>(static_cast<detid_t>(minval), static_cast<detid_t>(maxval),
                                                      time_conversion, m_calibration, m_scale_at_sample, m_masked);
    }
  }

  /// Copy histogram results and calculate errors
  void copyResultsAndCalculateErrors(API::ISpectrum &spectrum, const std::vector<uint32_t> &y_temp) const {
    auto &y_values = spectrum.dataY();
    std::copy(y_temp.cbegin(), y_temp.cend(), y_values.begin());
    auto &e_values = spectrum.dataE();
    std::transform(y_temp.cbegin(), y_temp.cend(), e_values.begin(),
                   [](uint32_t y) { return std::sqrt(static_cast<double>(y)); });
  }

  H5::H5File m_h5file;
  const std::vector<std::string> m_bankEntries;
  const std::map<detid_t, double> m_calibration;
  std::map<detid_t, double> m_scale_at_sample;
  const std::set<detid_t> m_masked;
  const size_t m_events_per_chunk;
  const size_t m_grainsize_event;
  std::shared_ptr<API::Progress> m_progress;

  // Logger reference (needs to be declared in implementation)
  static Kernel::Logger &g_log;
};

} // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
