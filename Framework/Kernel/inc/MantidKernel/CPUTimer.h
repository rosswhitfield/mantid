// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidKernel/DllConfig.h"
#include "MantidKernel/Timer.h"
#include <ctime>
#include <iosfwd>
#include <string>

namespace Mantid {
namespace Kernel {

/** CPUTimer : Timer that uses the CPU time, rather than wall-clock time
 * to measure execution time.
 *
 * @author Janik Zikovsky
 * @date 2011-04-04 12:17:48.579100
 */
class MANTID_KERNEL_DLL CPUTimer {
public:
  CPUTimer();
  float elapsedCPU(bool doReset = true);
  float elapsedWallClock(bool doReset = true);
  void reset();
  float CPUfraction(bool doReset = true);
  std::string str();

private:
  /// The starting time (implementation dependent format)
  clock_t m_start;

  /// The regular (wall-clock time).
  Timer m_wallClockTime;
};

MANTID_KERNEL_DLL std::ostream &operator<<(std::ostream &, CPUTimer &);

} // namespace Kernel
} // namespace Mantid
