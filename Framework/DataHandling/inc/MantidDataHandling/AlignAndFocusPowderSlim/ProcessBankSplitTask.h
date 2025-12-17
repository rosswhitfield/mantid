// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2025 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/MatrixWorkspace.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/NexusLoader.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/ProcessBankTaskBase.h"
#include "MantidGeometry/IDTypes.h"
#include <H5Cpp.h>
#include <MantidAPI/Progress.h>
#include <map>
#include <set>
#include <tbb/tbb.h>
#include <vector>

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

class ProcessBankSplitTask : public ProcessBankTaskBase {
public:
  ProcessBankSplitTask(std::vector<std::string> &bankEntryNames, H5::H5File &h5file, const bool is_time_filtered,
                       std::vector<int> &workspaceIndices, std::vector<API::MatrixWorkspace_sptr> &wksps,
                       const std::map<detid_t, double> &calibration, const std::map<detid_t, double> &scale_at_sample,
                       const std::set<detid_t> &masked, const size_t events_per_chunk, const size_t grainsize_event,
                       std::vector<std::pair<int, PulseROI>> target_to_pulse_indices,
                       std::shared_ptr<API::Progress> &progress);

  void operator()(const tbb::blocked_range<size_t> &range) const;

private:
  mutable NexusLoader m_loader;
  std::vector<int> m_workspaceIndices;
  std::vector<API::MatrixWorkspace_sptr> m_wksps;
};
} // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
