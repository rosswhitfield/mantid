// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#pragma once

#include "MantidDataHandling/DllConfig.h"
#include "MantidGeometry/IDTypes.h"
#include "MantidKernel/Matrix.h"

#include <vector>

namespace Mantid {
namespace Geometry {
class DetectorInfo;
}
namespace DataHandling::LoadEventNexusMD {

/**
 * Everything about a single detector that does not change from event to event. For elastic scattering with the beam
 * along the z-axis the momentum transfer of an event is
 *
 *   |k| = kfactor / tof     and     Q = |k| * (qx, qy, qz)
 *
 * with the time-of-flight in microseconds and |k| = 2*pi/lambda in inverse Angstrom. The direction is
 * qSign * (-ex, -ey, 1 - ez), where (ex, ey, ez) is the unit vector from the sample towards the detector, already
 * rotated into the target frame. Folding the rotation in here means Q_sample costs exactly as much per event as Q_lab.
 *
 * These are held in double rather than float so that the coordinates come out bit for bit identical to ConvertToMD,
 * which works in double throughout and rounds only when it stores the event.
 */
struct QCoefficients {
  double qx{0.};
  double qy{0.};
  double qz{0.};
  /// |k| in inverse Angstrom is this divided by the time-of-flight in microseconds. Zero marks an unusable detector.
  double kfactor{0.};
};

/// Returned for detector ids with no usable geometry: monitors, and gaps in the instrument's detector id range.
inline constexpr QCoefficients IGNORE_DETECTOR{};

/**
 * Conversion coefficients for a whole instrument, held in a vector offset by the smallest detector id because that is
 * a faster lookup than a map.
 *
 * The table is keyed by detector id rather than by bank because a NeXus bank name need not name an instrument
 * component at all: on CORELLI the panels are called A1 to C29 while the file still stores bank1_events. Keying on
 * the detector id is what LoadEventNexus does, and works whatever the instrument looks like.
 */
class MANTID_DATAHANDLING_DLL QConversionTable {
public:
  QConversionTable() = default;
  /**
   * @param detectorInfo Geometry for the whole instrument. Monitors are skipped.
   * @param rotation Applied to every direction. Identity for Q_lab, the inverse goniometer matrix for Q_sample.
   * @param qSign 1 for the "Inelastic" Q convention (ki-kf), -1 for "Crystallography".
   */
  QConversionTable(const Geometry::DetectorInfo &detectorInfo, const Kernel::DblMatrix &rotation, const double qSign);

  /// Coefficients for a detector id, or IGNORE_DETECTOR when no detector has that id
  const QCoefficients &value(const detid_t detid) const {
    if (detid < m_detid_offset)
      return IGNORE_DETECTOR;
    const auto index = static_cast<size_t>(detid - m_detid_offset);
    if (index >= m_coefficients.size())
      return IGNORE_DETECTOR;
    return m_coefficients[index];
  }

  bool empty() const { return m_coefficients.empty(); }
  detid_t idmin() const { return m_detid_offset; }
  detid_t idmax() const { return m_detid_offset + static_cast<detid_t>(m_coefficients.size()) - 1; }
  const std::vector<QCoefficients> &coefficients() const { return m_coefficients; }

private:
  std::vector<QCoefficients> m_coefficients;
  detid_t m_detid_offset{0};
};

} // namespace DataHandling::LoadEventNexusMD
} // namespace Mantid
