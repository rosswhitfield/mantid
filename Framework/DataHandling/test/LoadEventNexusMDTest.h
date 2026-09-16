// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include <cxxtest/TestSuite.h>

#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/ExperimentInfo.h"
#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/IMDEventWorkspace.h"
#include "MantidAPI/IMDHistoWorkspace.h"
#include "MantidAPI/Run.h"
#include "MantidDataHandling/LoadEventNexusMD.h"
#include "MantidGeometry/MDGeometry/QLab.h"
#include "MantidGeometry/MDGeometry/QSample.h"
#include "MantidKernel/SpecialCoordinateSystem.h"

#include <string>
#include <vector>

using Mantid::API::AnalysisDataService;
using Mantid::API::FrameworkManager;
using Mantid::API::IMDEventWorkspace_sptr;
using Mantid::API::IMDHistoWorkspace_sptr;
using Mantid::DataHandling::LoadEventNexusMD::LoadEventNexusMD;
using namespace Mantid::DataHandling::LoadEventNexusMD::PropertyNames;
namespace QFrames = Mantid::DataHandling::LoadEventNexusMD::QFrames;
namespace EventTypes = Mantid::DataHandling::LoadEventNexusMD::EventTypes;

namespace {
/// Small single crystal event file, 41 banks of a 256x256 MANDI detector
const std::string MANDI_FILE("MANDI_13064.nxs.h5");
/// a box that comfortably encloses every event in MANDI_FILE
const std::vector<double> MIN_VALS{-3., -3., -1.};
const std::vector<double> MAX_VALS{4., 5., 7.};
} // namespace

class LoadEventNexusMDTest : public CxxTest::TestSuite {
public:
  static LoadEventNexusMDTest *createSuite() { return new LoadEventNexusMDTest(); }
  static void destroySuite(LoadEventNexusMDTest *suite) { delete suite; }

  void test_Init() {
    LoadEventNexusMD alg;
    TS_ASSERT_THROWS_NOTHING(alg.initialize());
    TS_ASSERT(alg.isInitialized());
  }

  void test_extents_must_be_given_in_pairs() {
    assertRejected([](LoadEventNexusMD &alg) { alg.setProperty(MIN_VALUES, MIN_VALS); });
  }

  void test_extents_must_have_three_values() {
    assertRejected([](LoadEventNexusMD &alg) {
      alg.setProperty(MIN_VALUES, std::vector<double>{-1., -1.});
      alg.setProperty(MAX_VALUES, std::vector<double>{1., 1.});
    });
  }

  void test_extents_must_be_ordered() {
    assertRejected([](LoadEventNexusMD &alg) {
      alg.setProperty(MIN_VALUES, std::vector<double>{-1., -1., 2.});
      alg.setProperty(MAX_VALUES, std::vector<double>{1., 1., 1.});
    });
  }

  void test_chunk_must_not_be_smaller_than_the_grainsize() {
    assertRejected([](LoadEventNexusMD &alg) {
      alg.setProperty(READ_SIZE_FROM_DISK, 100);
      alg.setProperty(EVENTS_PER_THREAD, 1000);
    });
  }

  void test_exec_QSample() {
    const auto ws = runAlgorithm(QFrames::SAMPLE, EventTypes::FULL, MIN_VALS, MAX_VALS);

    TS_ASSERT_EQUALS(ws->getNumDims(), 3);
    TS_ASSERT_EQUALS(ws->getSpecialCoordinateSystem(), Mantid::Kernel::QSample);
    const std::vector<std::string> names{"Q_sample_x", "Q_sample_y", "Q_sample_z"};
    for (size_t d = 0; d < 3; ++d) {
      TS_ASSERT_EQUALS(ws->getDimension(d)->getName(), names[d]);
      TS_ASSERT_EQUALS(ws->getDimension(d)->getDimensionId(), "Q" + std::to_string(d + 1));
      TS_ASSERT_EQUALS(ws->getDimension(d)->getMDFrame().name(), Mantid::Geometry::QSample::QSampleName);
      TS_ASSERT_DELTA(ws->getDimension(d)->getMinimum(), MIN_VALS[d], 1e-5);
      TS_ASSERT_DELTA(ws->getDimension(d)->getMaximum(), MAX_VALS[d], 1e-5);
    }
    TS_ASSERT_EQUALS(ws->getNPoints(), 178721);

    // the instrument, logs and sample come across, and T0 is recorded as LoadEventNexus does
    TS_ASSERT_EQUALS(ws->getNumExperimentInfo(), 1);
    const auto &run = ws->getExperimentInfo(0)->run();
    TS_ASSERT(run.hasProperty("T0"));
    TS_ASSERT_DELTA(run.getPropertyAsSingleValue("T0"), -6.711, 1e-9);
    TS_ASSERT(run.hasProperty("W_MATRIX"));
  }

  void test_exec_QLab() {
    const auto ws = runAlgorithm(QFrames::LAB, EventTypes::FULL, MIN_VALS, MAX_VALS);
    TS_ASSERT_EQUALS(ws->getSpecialCoordinateSystem(), Mantid::Kernel::QLab);
    TS_ASSERT_EQUALS(ws->getDimension(0)->getName(), "Q_lab_x");
    TS_ASSERT_EQUALS(ws->getDimension(0)->getMDFrame().name(), Mantid::Geometry::QLab::QLabName);
    TS_ASSERT(ws->getNPoints() > 0);
  }

  /**
   * The load bearing test: for the same extents this must put every event in exactly the same place as
   * LoadEventNexus followed by ConvertToMD. It is what pins down the Q convention, the wavenumber constant and the
   * T0 correction.
   */
  void test_matches_LoadEventNexus_then_ConvertToMD_in_QSample() {
    assertMatchesConvertToMD(QFrames::SAMPLE, "Q_sample");
  }

  void test_matches_LoadEventNexus_then_ConvertToMD_in_QLab() { assertMatchesConvertToMD(QFrames::LAB, "Q_lab"); }

  void test_lean_events_hold_the_same_data() {
    const auto full = runAlgorithm(QFrames::SAMPLE, EventTypes::FULL, MIN_VALS, MAX_VALS, "md_full");
    const auto lean = runAlgorithm(QFrames::SAMPLE, EventTypes::LEAN, MIN_VALS, MAX_VALS, "md_lean");

    TS_ASSERT_EQUALS(full->id(), "MDEventWorkspace<MDEvent,3>");
    TS_ASSERT_EQUALS(lean->id(), "MDEventWorkspace<MDLeanEvent,3>");
    TS_ASSERT_EQUALS(lean->getNPoints(), full->getNPoints());
    assertSameBinnedSignal("md_full", "md_lean");

    AnalysisDataService::Instance().remove("md_full");
    AnalysisDataService::Instance().remove("md_lean");
  }

  void test_calculated_extents_keep_every_event() {
    // a deliberately over-wide box cannot reject anything, so the calculated one must find the same total
    const auto wide = runAlgorithm(QFrames::SAMPLE, EventTypes::FULL, std::vector<double>{-50., -50., -50.},
                                   std::vector<double>{50., 50., 50.});
    const auto calculated = runAlgorithm(QFrames::SAMPLE, EventTypes::FULL, {}, {});

    TS_ASSERT_EQUALS(calculated->getNPoints(), wide->getNPoints());
    // and the box must actually be tight around the data rather than the fallback
    for (size_t d = 0; d < 3; ++d) {
      TS_ASSERT(calculated->getDimension(d)->getMinimum() > -50.f);
      TS_ASSERT(calculated->getDimension(d)->getMaximum() < 50.f);
    }
  }

private:
  /// validateInputs is private, so bad input is checked through execute()
  template <typename Setup> void assertRejected(const Setup &setup) {
    LoadEventNexusMD alg;
    alg.setRethrows(true);
    alg.initialize();
    alg.setChild(true);
    alg.setPropertyValue(FILENAME, MANDI_FILE);
    alg.setPropertyValue(OUTPUT_WKSP, "unused_for_child");
    setup(alg);
    TS_ASSERT_THROWS(alg.execute(), const std::runtime_error &);
  }

  IMDEventWorkspace_sptr runAlgorithm(const std::string &qframe, const std::string &eventType,
                                      const std::vector<double> &minVals, const std::vector<double> &maxVals,
                                      const std::string &outName = "") {
    LoadEventNexusMD alg;
    alg.setRethrows(true);
    TS_ASSERT_THROWS_NOTHING(alg.initialize());
    alg.setChild(outName.empty());
    alg.setPropertyValue(FILENAME, MANDI_FILE);
    alg.setPropertyValue(OUTPUT_WKSP, outName.empty() ? "unused_for_child" : outName);
    alg.setPropertyValue(Q_FRAME, qframe);
    alg.setPropertyValue(EVENT_TYPE, eventType);
    if (!minVals.empty()) {
      alg.setProperty(MIN_VALUES, minVals);
      alg.setProperty(MAX_VALUES, maxVals);
    }
    TS_ASSERT_THROWS_NOTHING(alg.execute());
    TS_ASSERT(alg.isExecuted());
    // an output property releases its pointer once it has stored the workspace, so a named run comes back from the ADS
    if (!outName.empty())
      return AnalysisDataService::Instance().retrieveWS<Mantid::API::IMDEventWorkspace>(outName);
    return alg.getProperty(OUTPUT_WKSP);
  }

  /// Load the reference workspace once; it carries a 2.7 million detector instrument
  static std::string referenceEventWorkspace() {
    const std::string name("LoadEventNexusMDTest_events");
    if (!AnalysisDataService::Instance().doesExist(name))
      FrameworkManager::Instance().exec("LoadEventNexus", 4, "Filename", MANDI_FILE.c_str(), "OutputWorkspace",
                                        name.c_str());
    return name;
  }

  void assertMatchesConvertToMD(const std::string &qframe, const std::string &q3dFrame) {
    const std::string mine("LoadEventNexusMDTest_mine");
    const std::string reference("LoadEventNexusMDTest_reference");

    runAlgorithm(qframe, EventTypes::FULL, MIN_VALS, MAX_VALS, mine);

    FrameworkManager::Instance().exec("ConvertToMD", 14, "InputWorkspace", referenceEventWorkspace().c_str(),
                                      "OutputWorkspace", reference.c_str(), "QDimensions", "Q3D", "dEAnalysisMode",
                                      "Elastic", "Q3DFrames", q3dFrame.c_str(), "MinValues", "-3,-3,-1", "MaxValues",
                                      "4,5,7");

    const IMDEventWorkspace_sptr a = AnalysisDataService::Instance().retrieveWS<Mantid::API::IMDEventWorkspace>(mine);
    const IMDEventWorkspace_sptr b =
        AnalysisDataService::Instance().retrieveWS<Mantid::API::IMDEventWorkspace>(reference);
    TS_ASSERT_EQUALS(a->getNPoints(), b->getNPoints());
    assertSameBinnedSignal(mine, reference);

    AnalysisDataService::Instance().remove(mine);
    AnalysisDataService::Instance().remove(reference);
  }

  /// Bin both workspaces onto the same grid and require an identical count in every bin
  void assertSameBinnedSignal(const std::string &lhs, const std::string &rhs) {
    const auto binned = [](const std::string &in, const std::string &out) {
      FrameworkManager::Instance().exec("BinMD", 10, "InputWorkspace", in.c_str(), "OutputWorkspace", out.c_str(),
                                        "AlignedDim0", "Q1,-3,4,40", "AlignedDim1", "Q2,-3,5,40", "AlignedDim2",
                                        "Q3,-1,7,40");
      return AnalysisDataService::Instance().retrieveWS<Mantid::API::IMDHistoWorkspace>(out);
    };
    const auto ha = binned(lhs, "LoadEventNexusMDTest_binA");
    const auto hb = binned(rhs, "LoadEventNexusMDTest_binB");

    TS_ASSERT_EQUALS(ha->getNPoints(), hb->getNPoints());
    const auto *sa = ha->getSignalArray();
    const auto *sb = hb->getSignalArray();
    double worst{0.};
    for (size_t i = 0; i < ha->getNPoints(); ++i)
      worst = std::max(worst, std::abs(sa[i] - sb[i]));
    TSM_ASSERT_EQUALS("every bin must hold the same number of events", worst, 0.);

    AnalysisDataService::Instance().remove("LoadEventNexusMDTest_binA");
    AnalysisDataService::Instance().remove("LoadEventNexusMDTest_binB");
  }
};
