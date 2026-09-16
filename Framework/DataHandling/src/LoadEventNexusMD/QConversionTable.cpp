// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "MantidDataHandling/LoadEventNexusMD/QConversionTable.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"
#include "MantidKernel/PhysicalConstants.h"
#include "MantidKernel/V3D.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Mantid::DataHandling::LoadEventNexusMD {

namespace {
/**
 * Wavenumber in inverse Angstrom is this, times the total flight path in metres, divided by the time-of-flight in
 * microseconds. Using h_bar rather than h is what supplies the factor of 2*pi, so the result is 2*pi/lambda. Same
 * constant as ConvertToDiffractionMDWorkspace.
 */
constexpr double KFACTOR_PER_METRE = (PhysicalConstants::NeutronMass * 1e-10) / (1e-6 * PhysicalConstants::h_bar);
} // namespace

QConversionTable::QConversionTable(const Geometry::DetectorInfo &detectorInfo, const Kernel::DblMatrix &rotation,
                                   const double qSign) {
  // the coefficients are held in a vector offset by the smallest detector id, so find the range first
  detid_t detid_min{std::numeric_limits<detid_t>::max()};
  detid_t detid_max{std::numeric_limits<detid_t>::lowest()};
  for (size_t index = 0; index < detectorInfo.size(); ++index) {
    if (detectorInfo.isMonitor(index))
      continue;
    const auto detid = detectorInfo.detid(index);
    detid_min = std::min(detid_min, detid);
    detid_max = std::max(detid_max, detid);
  }
  if (detid_min > detid_max)
    return; // an instrument with nothing but monitors

  m_detid_offset = detid_min;
  m_coefficients.assign(static_cast<size_t>(detid_max - detid_min) + 1, IGNORE_DETECTOR);

  const double l1 = detectorInfo.l1();

  for (size_t index = 0; index < detectorInfo.size(); ++index) {
    if (detectorInfo.isMonitor(index))
      continue;

    const double twoTheta = detectorInfo.twoTheta(index);
    const double azimuthal = detectorInfo.azimuthal(index);
    const double sinTwoTheta = std::sin(twoTheta);

    // direction of the momentum transfer, k * (z_hat - e_hat) in the default "Inelastic" convention
    const Kernel::V3D direction =
        rotation * Kernel::V3D(qSign * -sinTwoTheta * std::cos(azimuthal), qSign * -sinTwoTheta * std::sin(azimuthal),
                               qSign * (1. - std::cos(twoTheta)));

    auto &coeff = m_coefficients[static_cast<size_t>(detectorInfo.detid(index) - m_detid_offset)];
    coeff.qx = direction.X();
    coeff.qy = direction.Y();
    coeff.qz = direction.Z();
    coeff.kfactor = KFACTOR_PER_METRE * (l1 + detectorInfo.l2(index));
  }
}

} // namespace Mantid::DataHandling::LoadEventNexusMD
