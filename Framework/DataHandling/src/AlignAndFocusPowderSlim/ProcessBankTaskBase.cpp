// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2025 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "MantidDataHandling/AlignAndFocusPowderSlim/ProcessBankTaskBase.h"
#include "MantidKernel/Logger.h"

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

namespace {
// Logger for base class
auto g_log_base = Kernel::Logger("ProcessBankTaskBase");
} // namespace

// Define static logger reference
Kernel::Logger &ProcessBankTaskBase::g_log = g_log_base;

} // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
