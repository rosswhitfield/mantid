// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include <cxxtest/TestSuite.h>

#include "MantidDataHandling/LoadEventNexusMD/QConversionTable.h"
#include "MantidFrameworkTestHelpers/WorkspaceCreationHelper.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"
#include "MantidKernel/Matrix.h"
#include "MantidKernel/UnitConversion.h"

#include <cmath>

using Mantid::detid_t;
using Mantid::DataHandling::LoadEventNexusMD::QConversionTable;
using Mantid::Kernel::DblMatrix;

class QConversionTableTest : public CxxTest::TestSuite {
public:
  static QConversionTableTest *createSuite() { return new QConversionTableTest(); }
  static void destroySuite(QConversionTableTest *suite) { delete suite; }

  /// 2 banks of 3x3 pixels, plus a monitor
  QConversionTableTest() : m_wksp(WorkspaceCreationHelper::create2DWorkspaceWithRectangularInstrument(2, 3, 1)) {}

  void test_covers_every_detector_id() {
    const QConversionTable table(detInfo(), identity(), 1.);
    TS_ASSERT(!table.empty());

    detid_t expectedMin{std::numeric_limits<detid_t>::max()};
    detid_t expectedMax{std::numeric_limits<detid_t>::lowest()};
    for (size_t i = 0; i < detInfo().size(); ++i) {
      if (detInfo().isMonitor(i))
        continue;
      expectedMin = std::min(expectedMin, detInfo().detid(i));
      expectedMax = std::max(expectedMax, detInfo().detid(i));
    }
    TS_ASSERT_EQUALS(table.idmin(), expectedMin);
    TS_ASSERT_EQUALS(table.idmax(), expectedMax);
  }

  void test_detids_outside_the_instrument_are_ignored() {
    const QConversionTable table(detInfo(), identity(), 1.);
    TS_ASSERT_EQUALS(table.value(table.idmin() - 1).kfactor, 0.);
    TS_ASSERT_EQUALS(table.value(table.idmax() + 1).kfactor, 0.);
    TS_ASSERT_EQUALS(table.value(1000000).kfactor, 0.);
  }

  /// The direction is qSign * (-ex, -ey, 1 - ez) for the unit vector (ex, ey, ez) sample -> detector
  void test_direction_matches_the_detector_geometry() {
    const QConversionTable table(detInfo(), identity(), 1.);

    size_t checked{0};
    for (size_t i = 0; i < detInfo().size(); ++i) {
      const auto detid = detInfo().detid(i);
      if (detInfo().isMonitor(i)) {
        TSM_ASSERT_EQUALS("monitors must be ignored", table.value(detid).kfactor, 0.);
        continue;
      }
      const double twoTheta = detInfo().twoTheta(i);
      const double azimuthal = detInfo().azimuthal(i);
      const auto &coeff = table.value(detid);
      TS_ASSERT_DELTA(coeff.qx, -std::sin(twoTheta) * std::cos(azimuthal), 1e-12);
      TS_ASSERT_DELTA(coeff.qy, -std::sin(twoTheta) * std::sin(azimuthal), 1e-12);
      TS_ASSERT_DELTA(coeff.qz, 1. - std::cos(twoTheta), 1e-12);
      ++checked;
    }
    TS_ASSERT(checked > 0);
  }

  /**
   * kfactor divided by the time-of-flight must give the same wavenumber as Mantid's own TOF to Momentum conversion.
   * This is what pins down the constant and the use of the total flight path.
   */
  void test_kfactor_matches_mantid_unit_conversion() {
    const QConversionTable table(detInfo(), identity(), 1.);
    const double tof{12345.};

    size_t checked{0};
    for (size_t i = 0; i < detInfo().size(); ++i) {
      if (detInfo().isMonitor(i))
        continue;
      const double expected = Mantid::Kernel::UnitConversion::run(
          "TOF", "Momentum", tof, detInfo().l1(), Mantid::Kernel::DeltaEMode::Elastic,
          {{Mantid::Kernel::UnitParams::l2, detInfo().l2(i)},
           {Mantid::Kernel::UnitParams::twoTheta, detInfo().twoTheta(i)}});
      TS_ASSERT_DELTA(table.value(detInfo().detid(i)).kfactor / tof, expected, 1e-9 * expected);
      ++checked;
    }
    TS_ASSERT(checked > 0);
  }

  void test_crystallography_convention_flips_the_direction() {
    const QConversionTable inelastic(detInfo(), identity(), 1.);
    const QConversionTable crystallography(detInfo(), identity(), -1.);

    for (detid_t detid = inelastic.idmin(); detid <= inelastic.idmax(); ++detid) {
      const auto &a = inelastic.value(detid);
      const auto &b = crystallography.value(detid);
      TS_ASSERT_DELTA(a.qx, -b.qx, 1e-12);
      TS_ASSERT_DELTA(a.qy, -b.qy, 1e-12);
      TS_ASSERT_DELTA(a.qz, -b.qz, 1e-12);
      TSM_ASSERT_EQUALS("the convention must not change the wavenumber", a.kfactor, b.kfactor);
    }
  }

  void test_rotation_is_applied_to_the_direction() {
    // quarter turn about z: (x, y, z) -> (-y, x, z)
    DblMatrix rotation(3, 3, false);
    rotation[0][1] = -1.;
    rotation[1][0] = 1.;
    rotation[2][2] = 1.;

    const QConversionTable unrotated(detInfo(), identity(), 1.);
    const QConversionTable rotated(detInfo(), rotation, 1.);

    for (detid_t detid = unrotated.idmin(); detid <= unrotated.idmax(); ++detid) {
      const auto &a = unrotated.value(detid);
      const auto &b = rotated.value(detid);
      TS_ASSERT_DELTA(b.qx, -a.qy, 1e-12);
      TS_ASSERT_DELTA(b.qy, a.qx, 1e-12);
      TS_ASSERT_DELTA(b.qz, a.qz, 1e-12);
      TS_ASSERT_EQUALS(a.kfactor, b.kfactor);
    }
  }

private:
  static DblMatrix identity() { return DblMatrix(3, 3, true); }
  const Mantid::Geometry::DetectorInfo &detInfo() const { return m_wksp->detectorInfo(); }

  Mantid::API::MatrixWorkspace_sptr m_wksp;
};
