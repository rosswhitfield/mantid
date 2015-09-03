#ifndef MANTID_ALGORITHMS_MAYERSSAMPLECORRECTIONTEST_H_
#define MANTID_ALGORITHMS_MAYERSSAMPLECORRECTIONTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidAlgorithms/SampleCorrections/MayersSampleCorrection.h"
#include <algorithm>
#include <cmath>

using Mantid::Algorithms::MayersSampleCorrection;

class MayersSampleCorrectionTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static MayersSampleCorrectionTest *createSuite() {
    return new MayersSampleCorrectionTest();
  }
  static void destroySuite(MayersSampleCorrectionTest *suite) { delete suite; }

  void test_attentuaton_correction_for_fixed_mur() {
    std::vector<double> dummy(2, 0.0);
    dummy[1] = 1.0;
    MayersSampleCorrection mscat(createTestParameters(), dummy, dummy, dummy);
    auto absFactor = mscat.calculateSelfAttenuation(0.01);

    const double delta = 1e-8;
    TS_ASSERT_DELTA(0.00030887, absFactor, delta);
  }

  // clang-format off
  void test_multiple_scattering_with_fixed_mur_and_absorption_correction_factor()
  // clang-format on
  {
    std::vector<double> dummy(2, 0.0);
    dummy[1] = 1.0;
    MayersSampleCorrection mscat(createTestParameters(), dummy, dummy, dummy);
    const size_t irp(1);
    const double muR(0.01), abs(0.0003);
    auto absFactor = mscat.calculateMS(irp, muR, abs);

    const double delta = 1e-8;
    TS_ASSERT_DELTA(0.00461391, absFactor.first, delta);
    TS_ASSERT_DELTA(67.25351289, absFactor.second, delta);
  }

  void test_corrects_both_absorption_and_multiple_scattering_for_point_data() {
    const size_t nypts(100);
    std::vector<double> signal(nypts, 2.0), tof(nypts), error(nypts);
    std::transform(signal.begin(), signal.end(), error.begin(), sqrt);
    double xcur(100.0);
    std::generate(tof.begin(), tof.end(), [&xcur] { return xcur++; });
    MayersSampleCorrection mscat(createTestParameters(), tof, signal, error);

    // Correct it
    mscat.apply(signal, error);

    // Check some values
    const double delta(1e-06);
    TS_ASSERT_DELTA(100.0, tof.front(), delta);
    TS_ASSERT_DELTA(199.0, tof.back(), delta);

    TS_ASSERT_DELTA(0.37497317, signal.front(), delta);
    TS_ASSERT_DELTA(0.37629282, signal.back(), delta);

    TS_ASSERT_DELTA(0.26514607, error.front(), delta);
    TS_ASSERT_DELTA(0.2660792, error.back(), delta);
  }

  // clang-format off
  void test_corrects_both_absorption_and_multiple_scattering_for_histogram_data()
  // clang-format on
  {
    const size_t nypts(100);
    std::vector<double> signal(nypts, 2.0), tof(nypts + 1), error(nypts);
    std::transform(signal.begin(), signal.end(), error.begin(), sqrt);
    // Generate a histogram with the same mid points as the point data example
    double xcur(99.5);
    std::generate(tof.begin(), tof.end(), [&xcur] {
      double xold = xcur;
      xcur += 1.0;
      return xold;
    });
    MayersSampleCorrection mscat(createTestParameters(), tof, signal, error);

    // Correct it
    mscat.apply(signal, error);

    // Check some values
    const double delta(1e-06);
    TS_ASSERT_DELTA(99.5, tof.front(), delta);
    TS_ASSERT_DELTA(199.5, tof.back(), delta);

    TS_ASSERT_DELTA(0.37497317, signal.front(), delta);
    TS_ASSERT_DELTA(0.37629282, signal.back(), delta);

    TS_ASSERT_DELTA(0.26514607, error.front(), delta);
    TS_ASSERT_DELTA(0.2660792, error.back(), delta);
  }

  // ---------------------- Failure tests -----------------------------
  void test_tof_not_monotonically_increasing_throws_invalid_argument() {
    const size_t nypts(10);
    std::vector<double> signal(nypts, 2.0), tof(nypts + 1), error(nypts);
    std::transform(signal.begin(), signal.end(), error.begin(), sqrt);
    double xcur(199.5);
    std::generate(tof.begin(), tof.end(), [&xcur] {
      double xold = xcur;
      xcur -= 1.0;
      return xold;
    });
    TS_ASSERT_THROWS(
        MayersSampleCorrection(createTestParameters(), tof, signal, error),
        std::invalid_argument);
  }

private:
  MayersSampleCorrection::Parameters createTestParameters() {
    // A bit like a POLARIS spectrum
    MayersSampleCorrection::Parameters pars;
    pars.l1 = 14.0;
    pars.l2 = 2.2;
    pars.twoTheta = 0.10821;
    pars.phi = 0.0;
    pars.rho = 0.07261;
    pars.sigmaSc = 5.1;
    pars.sigmaAbs = 5.08;
    pars.cylRadius = 0.0025;
    pars.cylHeight = 0.04;
    return pars;
  }
};

#endif /* MANTID_ALGORITHMS_MAYERSSAMPLECORRECTIONTEST_H_ */
