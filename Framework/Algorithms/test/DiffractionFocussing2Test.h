// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/FrameworkManager.h"
#include "MantidAPI/SpectraAxis.h"
#include "MantidAlgorithms/ConvertUnits.h"
#include "MantidAlgorithms/DiffractionFocussing2.h"
#include "MantidAlgorithms/MaskBins.h"
#include "MantidAlgorithms/Rebin.h"
#include "MantidDataHandling/ApplyDiffCal.h"
#include "MantidDataHandling/LoadNexus.h"
#include "MantidDataHandling/LoadRaw3.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidFrameworkTestHelpers/WorkspaceCreationHelper.h"
#include "MantidHistogramData/LinearGenerator.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidKernel/cow_ptr.h"
#include "Poco/StreamChannel.h"
#include <cxxtest/TestSuite.h>

using namespace Mantid;
using namespace Mantid::DataHandling;
using namespace Mantid::API;
using namespace Mantid::Kernel;
using namespace Mantid::Algorithms;
using namespace Mantid::DataObjects;
using Mantid::HistogramData::BinEdges;
using Mantid::Types::Event::TofEvent;

class DiffractionFocussing2Test : public CxxTest::TestSuite {
public:
  void testName() { TS_ASSERT_EQUALS(focus.name(), "DiffractionFocussing"); }

  void testVersion() { TS_ASSERT_EQUALS(focus.version(), 2); }

  void testInit() {
    focus.initialize();
    TS_ASSERT(focus.isInitialized());
  }

  void testExec() {
    Mantid::DataHandling::LoadNexus loader;
    loader.initialize();
    loader.setPropertyValue("Filename", "HRP38692a.nxs"); // HRP38692.raw spectrum range 320 to 330

    std::string outputSpace = "tofocus";
    loader.setPropertyValue("OutputWorkspace", outputSpace);
    TS_ASSERT_THROWS_NOTHING(loader.execute());
    TS_ASSERT(loader.isExecuted());

    Mantid::DataHandling::ApplyDiffCal applyDiffCal;
    applyDiffCal.setChild(true);
    applyDiffCal.initialize();
    applyDiffCal.setPropertyValue("InstrumentWorkspace", outputSpace);
    applyDiffCal.setPropertyValue("CalibrationFile", "hrpd_new_072_01.cal");
    TS_ASSERT_THROWS_NOTHING(applyDiffCal.execute());
    TS_ASSERT(applyDiffCal.isExecuted());

    // Have to convert because diffraction focussing wants d-spacing
    Mantid::Algorithms::ConvertUnits convert;
    convert.initialize();
    convert.setPropertyValue("InputWorkspace", outputSpace);
    convert.setPropertyValue("OutputWorkspace", outputSpace);
    convert.setPropertyValue("Target", "dSpacing");
    TS_ASSERT_THROWS_NOTHING(convert.execute());
    TS_ASSERT(convert.isExecuted());

    focus.setPropertyValue("InputWorkspace", outputSpace);
    focus.setPropertyValue("OutputWorkspace", "focusedWS");
    focus.setPropertyValue("GroupingFileName", "hrpd_new_072_01.cal");

    TS_ASSERT_THROWS_NOTHING(focus.execute());
    TS_ASSERT(focus.isExecuted());

    MatrixWorkspace_const_sptr output;
    TS_ASSERT_THROWS_NOTHING(output = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>("focusedWS"));

    // only 1 group for this limited range of spectra
    TS_ASSERT_EQUALS(output->getNumberHistograms(), 1);

    AnalysisDataService::Instance().remove(outputSpace);
    AnalysisDataService::Instance().remove("focusedWS");
  }

  void test_preserve_compare() {
    // processed nexus file in event mode
    const std::string PG3_FILE("PG3_46577.nxs.h5");

    const std::string inputWSname("PG3_46577");
    const std::string groupWS("PG3_group");
    { // load data from disk and get prepared for focussing
      auto loader = AlgorithmFactory::Instance().create("Load", -1); // "LoadNexusProcessed", -1);
      loader->initialize();
      loader->setProperty("Filename", PG3_FILE);
      loader->setProperty("OutputWorkspace", inputWSname);
      loader->execute();
      if (!loader->isExecuted())
        throw std::runtime_error("Failed to load " + PG3_FILE);

      auto cropworkspace = AlgorithmFactory::Instance().create("CropWorkspace", -1);
      cropworkspace->initialize();
      cropworkspace->setPropertyValue("InputWorkspace", inputWSname);
      cropworkspace->setPropertyValue("OutputWorkspace", inputWSname);
      cropworkspace->setProperty("XMin", 300.);
      cropworkspace->execute();
      if (!cropworkspace->isExecuted())
        throw std::runtime_error("Failed to CropWorkspace");

      auto compress = AlgorithmFactory::Instance().create("CompressEvents", -1);
      compress->initialize();
      compress->setPropertyValue("InputWorkspace", inputWSname);
      compress->setPropertyValue("OutputWorkspace", inputWSname);
      compress->execute();
      if (!compress->isExecuted())
        throw std::runtime_error("Failed to compress events");

      auto resamplex = AlgorithmFactory::Instance().create("ResampleX", -1);
      resamplex->initialize();
      resamplex->setPropertyValue("InputWorkspace", inputWSname);
      resamplex->setPropertyValue("OutputWorkspace", inputWSname);
      resamplex->setProperty("NumberBins", 6000);
      resamplex->setProperty("LogBinning", true);
      resamplex->execute();
      if (!resamplex->isExecuted())
        throw std::runtime_error("Failed to resample x-axis in TOF");

      auto convertunits = AlgorithmFactory::Instance().create("ConvertUnits", -1);
      convertunits->initialize();
      convertunits->setPropertyValue("InputWorkspace", inputWSname);
      convertunits->setPropertyValue("OutputWorkspace", inputWSname);
      convertunits->setProperty("Target", "dSpacing");
      convertunits->execute();
      if (!convertunits->isExecuted())
        throw std::runtime_error("Failed to convert to dspacing");

      auto creategroup = AlgorithmFactory::Instance().create("CreateGroupingWorkspace", -1);
      creategroup->initialize();
      creategroup->setPropertyValue("InputWorkspace", inputWSname);
      creategroup->setPropertyValue("OutputWorkspace", groupWS);
      creategroup->setProperty("GroupDetectorsBy", "All");
      creategroup->execute();
      if (!creategroup->isExecuted())
        throw std::runtime_error("Failed to create grouping workspace");
    }

    const std::string outputFalse(inputWSname + "_false");
    {
      DiffractionFocussing2 focussing;
      focussing.initialize();
      TS_ASSERT_THROWS_NOTHING(focussing.setPropertyValue("InputWorkspace", inputWSname));
      TS_ASSERT_THROWS_NOTHING(focussing.setPropertyValue("OutputWorkspace", outputFalse));
      TS_ASSERT_THROWS_NOTHING(focussing.setPropertyValue("GroupingWorkspace", groupWS));
      TS_ASSERT_THROWS_NOTHING(focussing.setProperty("PreserveEvents", false));
      TS_ASSERT_THROWS_NOTHING(focussing.execute());
      if (!focussing.isExecuted())
        throw std::runtime_error("Failed to DiffractionFocus PreserveEvents=False");
    }

    const std::string outputTrue(inputWSname + "_true");
    {
      DiffractionFocussing2 focussing;
      focussing.initialize();
      TS_ASSERT_THROWS_NOTHING(focussing.setPropertyValue("InputWorkspace", inputWSname));
      TS_ASSERT_THROWS_NOTHING(focussing.setPropertyValue("OutputWorkspace", outputTrue));
      TS_ASSERT_THROWS_NOTHING(focussing.setPropertyValue("GroupingWorkspace", groupWS));
      TS_ASSERT_THROWS_NOTHING(focussing.setProperty("PreserveEvents", true));
      TS_ASSERT_THROWS_NOTHING(focussing.execute());
      if (!focussing.isExecuted())
        throw std::runtime_error("Failed to DiffractionFocus PreserveEvents=True");
    }

    // cleanup inputs
    AnalysisDataService::Instance().remove(inputWSname);
    AnalysisDataService::Instance().remove(groupWS);

    // make sure workspace types do not match
    MatrixWorkspace_const_sptr wsFalse = AnalysisDataService::Instance().retrieveWS<const MatrixWorkspace>(outputFalse);
    MatrixWorkspace_const_sptr wsTrue = AnalysisDataService::Instance().retrieveWS<const MatrixWorkspace>(outputTrue);
    TS_ASSERT_DIFFERS(wsFalse->id(), wsTrue->id()); // different workspace types

    // check x-axis and total counts match - this is easier to debug than CompareWorkspaces
    const size_t NUM_HISTO = wsFalse->getNumberHistograms();
    for (size_t ws_index = 0; ws_index < NUM_HISTO; ++ws_index) {
      TS_ASSERT_EQUALS(wsFalse->readX(ws_index), wsTrue->readX(ws_index));
      const auto totalFalse = std::accumulate(wsFalse->readY(ws_index).cbegin(), wsFalse->readY(ws_index).cend(), 0.);
      const auto totalTrue = std::accumulate(wsTrue->readY(ws_index).cbegin(), wsTrue->readY(ws_index).cend(), 0.);
      TS_ASSERT_DELTA(totalFalse, totalTrue, .1);
    }

    // compare workspaces
    auto comparator = AlgorithmFactory::Instance().create("CompareWorkspaces", -1);
    comparator->initialize();
    comparator->setProperty("Workspace1", outputFalse);
    comparator->setProperty("Workspace2", outputTrue);
    comparator->setProperty("CheckType", false);
    comparator->execute();
    bool result = comparator->getProperty("Result");
    TS_ASSERT(result);

    // cleanup outputs
    AnalysisDataService::Instance().remove(outputFalse);
    AnalysisDataService::Instance().remove(outputTrue);
  }

  void test_EventWorkspace_SameOutputWS() { dotestEventWorkspace(true, 2); }

  void test_EventWorkspace_DifferentOutputWS() { dotestEventWorkspace(false, 2); }
  void test_EventWorkspace_SameOutputWS_oneGroup() { dotestEventWorkspace(true, 1); }

  void test_EventWorkspace_DifferentOutputWS_oneGroup() { dotestEventWorkspace(false, 1); }

  void test_EventWorkspace_TwoGroups_dontPreserveEvents() { dotestEventWorkspace(false, 2, false); }

  void test_EventWorkspace_OneGroup_dontPreserveEvents() { dotestEventWorkspace(false, 1, false); }

  void test_rebin_parameters_histogram() {
    std::string inputWS("DiffractionFocussing2TestParam_ws");
    std::string groupWS("DiffractionFocussing2TestParam_groups");

    // histogram input
    const double xMin = 200.0;
    const double xMax = 600.0;
    create_test_workspace_input(inputWS, groupWS, xMin, xMax, "Histogram");
    run_rebin_parameters_test(inputWS, groupWS, "Workspace2D", xMin, xMax, true);
  }

  void test_rebin_parameters_event_preserve_events() {
    std::string inputWS("DiffractionFocussing2TestParam_ws");
    std::string groupWS("DiffractionFocussing2TestParam_groups");

    // histogram input
    // event input, PreserveEvent=True
    const double xMin = 200.0;
    const double xMax = 600.0;
    create_test_workspace_input(inputWS, groupWS, xMin, xMax, "Event");
    run_rebin_parameters_test(inputWS, groupWS, "EventWorkspace", xMin, xMax, true);
  }

  void test_rebin_parameters_event_dont_preserve_events() {
    std::string inputWS("DiffractionFocussing2TestParam_ws");
    std::string groupWS("DiffractionFocussing2TestParam_groups");

    // event input, PreserveEvent=False
    const double xMin = 200.0;
    const double xMax = 600.0;
    create_test_workspace_input(inputWS, groupWS, xMin, xMax, "Event");
    run_rebin_parameters_test(inputWS, groupWS, "Workspace2D", xMin, xMax, false);
  }

  void create_test_workspace_input(std::string inputWS, std::string groupWS, double xMin, double xMax,
                                   std::string workspaceType) {
    auto createWS = AlgorithmFactory::Instance().create("CreateSampleWorkspace", -1);
    createWS->initialize();
    createWS->setProperty("WorkspaceType", workspaceType);
    createWS->setProperty("XUnit", "dSpacing");
    createWS->setProperty("NumBanks", 2);
    createWS->setProperty("BankPixelWidth", 2);
    createWS->setProperty("Function", "Flat background");
    createWS->setProperty("XMin", xMin);
    createWS->setProperty("XMax", xMax);
    createWS->setProperty("BinWidth", 10.);
    createWS->setProperty("NumEvents", 40);
    createWS->setPropertyValue("OutputWorkspace", inputWS);
    createWS->execute();

    auto createGroupWS = AlgorithmFactory::Instance().create("CreateGroupingWorkspace", -1);
    createGroupWS->initialize();
    createGroupWS->setPropertyValue("InputWorkspace", inputWS);
    createGroupWS->setPropertyValue("GroupDetectorsBy", "bank");
    createGroupWS->setPropertyValue("OutputWorkspace", groupWS);
    createGroupWS->execute();
  }

  void run_rebin_parameters_test(std::string inputWS, std::string groupWS, std::string outputType, double xMin,
                                 double xMax, bool preserveEvents) {

    {
      const double test_xMin = xMin;
      const double test_xMax = xMax - 200.0;
      const double test_delta = 100.0;

      std::vector<double> dmins = {test_xMin};
      std::vector<double> dmaxs = {test_xMax};
      std::vector<double> delta = {test_delta};
      auto output = run_DiffractionFocussing2(inputWS, groupWS, dmins, dmaxs, delta, preserveEvents);
      if (!output)
        return;

      TS_ASSERT_EQUALS(output->id(), outputType)

      TS_ASSERT_EQUALS(output->getNumberHistograms(), 2);
      auto X0 = output->x(0);
      TS_ASSERT_EQUALS(X0.size(), 3.);
      TS_ASSERT_EQUALS(X0.front(), test_xMin);
      TS_ASSERT_EQUALS(X0.back(), test_xMax);
      auto X1 = output->x(1);
      TS_ASSERT_EQUALS(X1.size(), 3.);
      TS_ASSERT_EQUALS(X1.front(), test_xMin);
      TS_ASSERT_EQUALS(X1.back(), test_xMax);

      auto Y0 = output->y(0);
      TS_ASSERT_EQUALS(Y0.size(), 2.);
      TS_ASSERT_EQUALS(Y0[0], 40.);
      TS_ASSERT_EQUALS(Y0[1], 40.);
      auto Y1 = output->y(1);
      TS_ASSERT_EQUALS(Y1.size(), 2.);
      TS_ASSERT_EQUALS(Y1[0], 40.);
      TS_ASSERT_EQUALS(Y1[1], 40.);
    }

    {
      const double test_xMin = xMin;
      const double test_xMax = xMax - 200.0;
      const double test_delta = 100.0;

      std::vector<double> dmins = {test_xMin, test_xMin + 100.0};
      std::vector<double> dmaxs = {test_xMax};
      std::vector<double> delta = {test_delta};
      auto output = run_DiffractionFocussing2(inputWS, groupWS, dmins, dmaxs, delta, preserveEvents);
      if (!output)
        return;

      TS_ASSERT_EQUALS(output->id(), outputType)

      TS_ASSERT_EQUALS(output->getNumberHistograms(), 2);
      auto X0 = output->x(0);
      TS_ASSERT_EQUALS(X0.size(), 3.);
      TS_ASSERT_EQUALS(X0.front(), test_xMin);
      TS_ASSERT_EQUALS(X0.back(), test_xMax);
      auto X1 = output->x(1);
      TS_ASSERT_EQUALS(X1.size(), 2.);
      TS_ASSERT_EQUALS(X1.front(), test_xMin + 100.0);
      TS_ASSERT_EQUALS(X1.back(), test_xMax);
    }

    {
      const double test_xMin = xMin;
      const double test_xMax = xMax - 200.0;
      const double test_delta = 100.0;

      std::vector<double> dmins = {test_xMin};
      std::vector<double> dmaxs = {test_xMax};
      std::vector<double> delta = {test_delta, test_delta + 100.0};
      auto output = run_DiffractionFocussing2(inputWS, groupWS, dmins, dmaxs, delta, preserveEvents);
      if (!output)
        return;

      TS_ASSERT_EQUALS(output->id(), outputType)

      TS_ASSERT_EQUALS(output->getNumberHistograms(), 2);
      auto X0 = output->x(0);
      TS_ASSERT_EQUALS(X0.size(), 3.);
      TS_ASSERT_EQUALS(X0.front(), test_xMin);
      TS_ASSERT_EQUALS(X0.back(), test_xMax);
      auto X1 = output->x(1);
      TS_ASSERT_EQUALS(X1.size(), 2.);
      TS_ASSERT_EQUALS(X1.front(), test_xMin);
      TS_ASSERT_EQUALS(X1.back(), test_xMax);
    }

    {
      const double test_xMin = xMin;
      const double test_xMax = xMax;
      const double test_delta = 100.0;

      std::vector<double> dmins = {test_xMin + 100.0, test_xMin};
      std::vector<double> dmaxs = {test_xMax - 100.0, test_xMax};
      std::vector<double> delta = {test_delta, test_delta + 100.0};
      auto output = run_DiffractionFocussing2(inputWS, groupWS, dmins, dmaxs, delta, preserveEvents);
      if (!output)
        return;

      TS_ASSERT_EQUALS(output->id(), outputType)

      TS_ASSERT_EQUALS(output->getNumberHistograms(), 2);
      auto X0 = output->x(0);
      TS_ASSERT_EQUALS(X0.size(), 3.);
      TS_ASSERT_EQUALS(X0.front(), test_xMin + 100.0);
      TS_ASSERT_EQUALS(X0.back(), test_xMax - 100.0);
      auto X1 = output->x(1);
      TS_ASSERT_EQUALS(X1.size(), 3.);
      TS_ASSERT_EQUALS(X1.front(), test_xMin);
      TS_ASSERT_EQUALS(X1.back(), test_xMax);
    }

    {
      // Verify case where input data domain is contained by output domain:
      //   ensure that there are no NaN in the output workspace.  This tests a bug fix to
      //   an issue caused by the fact that in this case, the rebinning weights near the boundaries
      //   were previously set to zero.
      //  (See: EWM defect #12058 [ORNL])

      const double test_xMin = xMin;
      const double test_xMax = xMax;
      const double test_delta = 50.0;

      std::vector<double> dmins = {test_xMin - 100.0, test_xMin};
      std::vector<double> dmaxs = {test_xMax + 100.0, test_xMax};
      std::vector<double> delta = {test_delta, test_delta};
      auto output = run_DiffractionFocussing2(inputWS, groupWS, dmins, dmaxs, delta, preserveEvents);
      if (!output)
        return;

      TS_ASSERT_EQUALS(output->id(), outputType)

      TS_ASSERT_EQUALS(output->getNumberHistograms(), 2);
      const auto &X0 = output->x(0);
      const auto &Y0 = output->y(0);
      TS_ASSERT_EQUALS(X0.size(), 13);
      TS_ASSERT_EQUALS(X0.front(), test_xMin - 100.0);
      TS_ASSERT_EQUALS(X0.back(), test_xMax + 100.0);

      size_t nan_count = 0;
      for (const auto &y : Y0)
        if (std::isnan(y))
          ++nan_count;
      TS_ASSERT_EQUALS(nan_count, 0);

      const auto &X1 = output->x(1);
      TS_ASSERT_EQUALS(X1.size(), 9);
      TS_ASSERT_EQUALS(X1.front(), test_xMin);
      TS_ASSERT_EQUALS(X1.back(), test_xMax);
    }

    {
      // Test `FullBinsOnly=true` case.
      // (Previous test cases check `FullBinsOnly=false`, which is the default setting.)

      const double test_xMin = xMin;
      const double test_xMax = xMax;
      const double test_delta = 100.0;

      std::vector<double> dmins = {test_xMin + 100.0, test_xMin};
      std::vector<double> dmaxs = {test_xMax - 50.0, test_xMax};
      std::vector<double> delta = {test_delta, test_delta + 100.0};
      auto output = run_DiffractionFocussing2(inputWS, groupWS, dmins, dmaxs, delta, preserveEvents, true);
      if (!output)
        return;

      TS_ASSERT_EQUALS(output->id(), outputType)

      TS_ASSERT_EQUALS(output->getNumberHistograms(), 2);

      // First spectrum should not have its final bin.
      auto X0 = output->x(0);
      TS_ASSERT_EQUALS(X0.size(), 3.);
      TS_ASSERT_EQUALS(X0.front(), test_xMin + 100.0);
      TS_ASSERT_EQUALS(X0.back(), test_xMax - 100.0);

      // Second spectrum should be as before.
      auto X1 = output->x(1);
      TS_ASSERT_EQUALS(X1.size(), 3.);
      TS_ASSERT_EQUALS(X1.front(), test_xMin);
      TS_ASSERT_EQUALS(X1.back(), test_xMax);
    }

    {
      // Test `FullBinsOnly` case: negative test.

      const double test_xMin = xMin;
      const double test_xMax = xMax;
      const double test_delta = 100.0;

      std::vector<double> dmins = {test_xMin + 100.0, test_xMin};
      std::vector<double> dmaxs = {test_xMax - 50.0, test_xMax};
      std::vector<double> delta = {test_delta, test_delta + 100.0};
      auto output = run_DiffractionFocussing2(inputWS, groupWS, dmins, dmaxs, delta, preserveEvents);
      if (!output)
        return;

      TS_ASSERT_EQUALS(output->id(), outputType)

      TS_ASSERT_EQUALS(output->getNumberHistograms(), 2);

      // First spectrum should retain its final, partial bin.
      auto X0 = output->x(0);
      TS_ASSERT_EQUALS(X0.size(), 4.);
      TS_ASSERT_EQUALS(X0.front(), test_xMin + 100.0);
      TS_ASSERT_EQUALS(X0.back(), test_xMax - 50.0);

      // Second spectrum should be as before.
      auto X1 = output->x(1);
      TS_ASSERT_EQUALS(X1.size(), 3.);
      TS_ASSERT_EQUALS(X1.front(), test_xMin);
      TS_ASSERT_EQUALS(X1.back(), test_xMax);
    }

    // failure cases

    // NaN or delta=0
    run_DiffractionFocussing2(inputWS, groupWS, {std::numeric_limits<double>::quiet_NaN()}, {400}, {100},
                              preserveEvents, false, true);
    run_DiffractionFocussing2(inputWS, groupWS, {200}, {std::numeric_limits<double>::quiet_NaN()}, {100},
                              preserveEvents, false, true);
    run_DiffractionFocussing2(inputWS, groupWS, {200}, {400}, {std::numeric_limits<double>::quiet_NaN()},
                              preserveEvents, false, true);
    run_DiffractionFocussing2(inputWS, groupWS, {200}, {400}, {0}, preserveEvents, false, true);

    // must define all or none binning params
    run_DiffractionFocussing2(inputWS, groupWS, {}, {400}, {100}, preserveEvents, false, true);
    run_DiffractionFocussing2(inputWS, groupWS, {200}, {}, {100}, preserveEvents, false, true);
    run_DiffractionFocussing2(inputWS, groupWS, {200}, {400}, {}, preserveEvents, false, true);

    // dmax is not larger than dmin
    run_DiffractionFocussing2(inputWS, groupWS, {200, 400}, {400}, {100}, preserveEvents, false, true);
    run_DiffractionFocussing2(inputWS, groupWS, {200}, {200, 400}, {100}, preserveEvents, false, true);
    run_DiffractionFocussing2(inputWS, groupWS, {200, 400}, {300, 400}, {100}, preserveEvents, false, true);

    // too many parameters
    run_DiffractionFocussing2(inputWS, groupWS, {200, 200, 200}, {400}, {100}, preserveEvents, false, true);
    run_DiffractionFocussing2(inputWS, groupWS, {200}, {400, 400, 400}, {100}, preserveEvents, false, true);
    run_DiffractionFocussing2(inputWS, groupWS, {200}, {400}, {100, 100, 100}, preserveEvents, false, true);

    // cleanup workspaces
    AnalysisDataService::Instance().remove(inputWS);
    AnalysisDataService::Instance().remove(groupWS);
  }

  MatrixWorkspace_const_sptr run_DiffractionFocussing2(std::string inputWS, std::string groupWS,
                                                       const std::vector<double> dmins, const std::vector<double> dmaxs,
                                                       const std::vector<double> delta, bool preserveEvents,
                                                       bool fullBinsOnly = false, bool expectFailure = false) {
    std::string outputWS = inputWS + "_focussed";
    MatrixWorkspace_const_sptr output;

    DiffractionFocussing2 focus;
    focus.initialize();
    TS_ASSERT_THROWS_NOTHING(focus.setPropertyValue("InputWorkspace", inputWS));
    TS_ASSERT_THROWS_NOTHING(focus.setPropertyValue("OutputWorkspace", outputWS));
    TS_ASSERT_THROWS_NOTHING(focus.setPropertyValue("GroupingWorkspace", groupWS));
    TS_ASSERT_THROWS_NOTHING(focus.setProperty("PreserveEvents", preserveEvents));
    TS_ASSERT_THROWS_NOTHING(focus.setProperty("DMin", dmins));
    TS_ASSERT_THROWS_NOTHING(focus.setProperty("DMax", dmaxs));
    TS_ASSERT_THROWS_NOTHING(focus.setProperty("Delta", delta));
    TS_ASSERT_THROWS_NOTHING(focus.setProperty("FullBinsOnly", fullBinsOnly));

    if (expectFailure) {
      try {
        focus.execute();
      } catch (std::runtime_error &) {
      }
      TS_ASSERT(!focus.isExecuted());
    } else {
      TS_ASSERT_THROWS_NOTHING(focus.execute(););
      TS_ASSERT(focus.isExecuted());
      TS_ASSERT_THROWS_NOTHING(output = AnalysisDataService::Instance().retrieveWS<const MatrixWorkspace>(outputWS));
      TS_ASSERT(output);
    }
    return output;
  }

  void test_tof_deprecation_error_thrown() {
    // Simple test to be removed when TOF support is removed
    std::string nxsWSname("DiffractionFocussing2Test_ws");
    // Create the fake event workspace
    EventWorkspace_sptr inputW = WorkspaceCreationHelper::createEventWorkspaceWithFullInstrument(3, 1);
    AnalysisDataService::Instance().addOrReplace(nxsWSname, inputW);
    // Set xunit TOF
    inputW->getAxis(0)->unit() = UnitFactory::Instance().create("TOF");
    // Create a grouping workspace
    std::string GroupNames = "bank3";
    std::string groupWSName("DiffractionFocussing2Test_group");
    FrameworkManager::Instance().exec("CreateGroupingWorkspace", 6, "InputWorkspace", nxsWSname.c_str(), "GroupNames",
                                      GroupNames.c_str(), "OutputWorkspace", groupWSName.c_str());
    // Run algorithm
    DiffractionFocussing2 focus;
    focus.initialize();
    TS_ASSERT_THROWS_NOTHING(focus.setPropertyValue("InputWorkspace", nxsWSname));
    TS_ASSERT_THROWS_NOTHING(focus.setPropertyValue("OutputWorkspace", nxsWSname));
    TS_ASSERT_THROWS_NOTHING(focus.setPropertyValue("GroupingWorkspace", groupWSName));
    TS_ASSERT_THROWS_NOTHING(focus.setProperty("PreserveEvents", false));
    // setup poco stream to catch log output
    std::ostringstream oss;
    auto psc = new Poco::StreamChannel(oss);
    Poco::Logger::setChannel("DiffractionFocussing", psc);
    TS_ASSERT_THROWS_NOTHING(focus.execute(););
    // assert output contains deprecation message
    const auto logMsg = oss.str();
    TS_ASSERT(logMsg.find("Support for TOF data in DiffractionFocussing is deprecated") != std::string::npos)
  }

  void dotestEventWorkspace(bool inplace, size_t numgroups, bool preserveEvents = true, int bankWidthInPixels = 16) {
    std::string nxsWSname("DiffractionFocussing2Test_ws");

    // Create the fake event workspace
    EventWorkspace_sptr inputW = WorkspaceCreationHelper::createEventWorkspaceWithFullInstrument(3, bankWidthInPixels);
    AnalysisDataService::Instance().addOrReplace(nxsWSname, inputW);

    //    //----- Load some event data --------
    //    FrameworkManager::Instance().exec("LoadEventNexus", 4,
    //        "Filename", "CNCS_7860_event.nxs",
    //        "OutputWorkspace", nxsWSname.c_str());

    //-------- Check on the input workspace ---------------
    TS_ASSERT(inputW);
    if (!inputW)
      return;

    // Fake a d-spacing unit in the data.
    inputW->getAxis(0)->unit() = UnitFactory::Instance().create("dSpacing");

    // Create a DIFFERENT x-axis for each pixel. Starting bin = the input
    // workspace index #
    for (size_t pix = 0; pix < inputW->getNumberHistograms(); pix++) {
      // Set an X-axis
      double x = static_cast<double>(1 + pix);
      inputW->setHistogram(pix, BinEdges{x + 0, x + 1, x + 2, x + 3, 1e6});
      inputW->getSpectrum(pix).addEventQuickly(TofEvent(1000.0));
    }

    // ------------ Create a grouping workspace by name -------------
    std::string GroupNames = "bank2,bank3";
    if (numgroups == 1)
      GroupNames = "bank3";
    std::string groupWSName("DiffractionFocussing2Test_group");
    FrameworkManager::Instance().exec("CreateGroupingWorkspace", 6, "InputWorkspace", nxsWSname.c_str(), "GroupNames",
                                      GroupNames.c_str(), "OutputWorkspace", groupWSName.c_str());

    // ------------ Create a grouping workspace by name -------------
    DiffractionFocussing2 focus;
    focus.initialize();
    TS_ASSERT_THROWS_NOTHING(focus.setPropertyValue("InputWorkspace", nxsWSname));
    std::string outputws = nxsWSname + "_focussed";
    if (inplace)
      outputws = nxsWSname;
    TS_ASSERT_THROWS_NOTHING(focus.setPropertyValue("OutputWorkspace", outputws));

    // This fake calibration file was generated using
    // DiffractionFocussing2Test_helper.py
    TS_ASSERT_THROWS_NOTHING(focus.setPropertyValue("GroupingWorkspace", groupWSName));
    TS_ASSERT_THROWS_NOTHING(focus.setProperty("PreserveEvents", preserveEvents));
    // OK, run the algorithm
    TS_ASSERT_THROWS_NOTHING(focus.execute(););
    TS_ASSERT(focus.isExecuted());

    MatrixWorkspace_const_sptr output;
    TS_ASSERT_THROWS_NOTHING(output = AnalysisDataService::Instance().retrieveWS<const MatrixWorkspace>(outputws));
    if (!output)
      return;

    // ---- Did we keep the event workspace ----
    EventWorkspace_const_sptr outputEvent;
    TS_ASSERT_THROWS_NOTHING(outputEvent = std::dynamic_pointer_cast<const EventWorkspace>(output));
    if (preserveEvents) {
      TS_ASSERT(outputEvent);
      if (!outputEvent)
        return;
    } else {
      TS_ASSERT(!outputEvent);
    }

    TS_ASSERT_EQUALS(output->getNumberHistograms(), numgroups);
    if (output->getNumberHistograms() != numgroups)
      return;

    TS_ASSERT_EQUALS(output->blocksize(), 4);

    TS_ASSERT_EQUALS(output->getAxis(1)->length(), numgroups);
    TS_ASSERT_EQUALS(output->getSpectrum(0).getSpectrumNo(), 1);

    // Events in these two banks alone
    if (preserveEvents)
      TS_ASSERT_EQUALS(outputEvent->getNumberEvents(), (numgroups == 2) ? (bankWidthInPixels * bankWidthInPixels * 2)
                                                                        : bankWidthInPixels * bankWidthInPixels);

    // Now let's test the grouping of detector UDETS to groups
    for (size_t wi = 0; wi < output->getNumberHistograms(); wi++) {
      // This is the list of the detectors (grouped)
      auto mylist = output->getSpectrum(wi).getDetectorIDs();
      // 1024 pixels in a bank
      TS_ASSERT_EQUALS(mylist.size(), bankWidthInPixels * bankWidthInPixels);
    }

    if (preserveEvents) {
      // Now let's try to rebin using log parameters (this used to fail?)
      Rebin rebin;
      rebin.initialize();
      rebin.setPropertyValue("InputWorkspace", outputws);
      rebin.setPropertyValue("OutputWorkspace", outputws);
      // Check it fails if "Params" property not set
      rebin.setPropertyValue("Params", "2.0,-1.0,65535");
      TS_ASSERT(rebin.execute());
      TS_ASSERT(rebin.isExecuted());

      /* Get the output ws again */
      outputEvent = AnalysisDataService::Instance().retrieveWS<EventWorkspace>(outputws);
      double events_after_binning = 0;
      for (size_t workspace_index = 0; workspace_index < outputEvent->getNumberHistograms(); workspace_index++) {
        // should be 16 bins
        TS_ASSERT_EQUALS(outputEvent->refX(workspace_index)->size(), 16);
        // There should be some data in the bins
        for (int i = 0; i < 15; i++)
          events_after_binning += outputEvent->dataY(workspace_index)[i];
      }
      // The count sums up to the same as the number of events
      TS_ASSERT_DELTA(events_after_binning,
                      (numgroups == 2) ? double(bankWidthInPixels * bankWidthInPixels) * 2.0
                                       : double(bankWidthInPixels * bankWidthInPixels),
                      1e-4);
    }
  }

private:
  DiffractionFocussing2 focus;
};

//================================================================================================
//================================================================================================
//================================================================================================

class DiffractionFocussing2TestPerformance : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static DiffractionFocussing2TestPerformance *createSuite() { return new DiffractionFocussing2TestPerformance(); }
  static void destroySuite(DiffractionFocussing2TestPerformance *suite) { delete suite; }

  EventWorkspace_sptr ws;

  DiffractionFocussing2TestPerformance() {
    auto alg = AlgorithmFactory::Instance().create("LoadEmptyInstrument", 1);
    alg->initialize();
    alg->setRethrows(true);
    alg->setPropertyValue("Filename", "SNAP_Definition_2011-09-07.xml");
    alg->setPropertyValue("OutputWorkspace", "SNAP_empty");
    alg->setPropertyValue("MakeEventWorkspace", "1");
    alg->execute();
    ws = AnalysisDataService::Instance().retrieveWS<EventWorkspace>("SNAP_empty");
    ws->sortAll(TOF_SORT, nullptr);

    // Fill a whole bunch of events
    PARALLEL_FOR_NO_WSP_CHECK()
    for (int i = 0; i < static_cast<int>(ws->getNumberHistograms()); i++) {
      EventList &el = ws->getSpectrum(i);
      for (int j = 0; j < 20; j++) {
        el.addEventQuickly(TofEvent(double(j) * 1e-3));
      }
    }
    ws->getAxis(0)->setUnit("dSpacing");
    // Create the x-axis for histogramming.
    ws->setAllX(BinEdges{0.0, 1e6});

    alg = AlgorithmFactory::Instance().create("CreateGroupingWorkspace", 1);
    alg->initialize();
    alg->setPropertyValue("InputWorkspace", "SNAP_empty");
    alg->setPropertyValue("GroupNames", "bank1");
    alg->setPropertyValue("OutputWorkspace", "SNAP_group_bank1");
    alg->execute();

    alg = AlgorithmFactory::Instance().create("CreateGroupingWorkspace", 1);
    alg->initialize();
    alg->setPropertyValue("InputWorkspace", "SNAP_empty");
    alg->setPropertyValue("GroupNames", "bank1,bank2,bank3,bank4,bank5,bank6");
    alg->setPropertyValue("OutputWorkspace", "SNAP_group_several");
    alg->execute();
  }

  ~DiffractionFocussing2TestPerformance() override {
    AnalysisDataService::Instance().remove("SNAP_empty");
    AnalysisDataService::Instance().remove("SNAP_group_bank1");
    AnalysisDataService::Instance().remove("SNAP_group_several");
  }

  void test_SNAP_event_one_group() {
    auto alg = AlgorithmFactory::Instance().create("DiffractionFocussing", 2);
    alg->initialize();
    alg->setPropertyValue("InputWorkspace", "SNAP_empty");
    alg->setPropertyValue("GroupingWorkspace", "SNAP_group_bank1");
    alg->setPropertyValue("OutputWorkspace", "SNAP_focus");
    alg->setPropertyValue("PreserveEvents", "1");
    alg->execute();
    EventWorkspace_sptr outWS = AnalysisDataService::Instance().retrieveWS<EventWorkspace>("SNAP_focus");

    TS_ASSERT_EQUALS(outWS->getNumberHistograms(), 1);
    TS_ASSERT_EQUALS(outWS->getNumberEvents(), 20 * 65536);
    AnalysisDataService::Instance().remove("SNAP_focus");
  }

  void test_SNAP_event_six_groups() {
    auto alg = AlgorithmFactory::Instance().create("DiffractionFocussing", 2);
    alg->initialize();
    alg->setPropertyValue("InputWorkspace", "SNAP_empty");
    alg->setPropertyValue("GroupingWorkspace", "SNAP_group_several");
    alg->setPropertyValue("OutputWorkspace", "SNAP_focus");
    alg->setPropertyValue("PreserveEvents", "1");
    alg->execute();
    EventWorkspace_sptr outWS = AnalysisDataService::Instance().retrieveWS<EventWorkspace>("SNAP_focus");

    TS_ASSERT_EQUALS(outWS->getNumberHistograms(), 6);
    TS_ASSERT_EQUALS(outWS->getNumberEvents(), 6 * 20 * 65536);
    AnalysisDataService::Instance().remove("SNAP_focus");
  }

  void test_SNAP_event_one_group_dontPreserveEvents() {
    auto alg = AlgorithmFactory::Instance().create("DiffractionFocussing", 2);
    alg->initialize();
    alg->setPropertyValue("InputWorkspace", "SNAP_empty");
    alg->setPropertyValue("GroupingWorkspace", "SNAP_group_bank1");
    alg->setPropertyValue("OutputWorkspace", "SNAP_focus");
    alg->setPropertyValue("PreserveEvents", "0");
    alg->execute();
    MatrixWorkspace_sptr outWS = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>("SNAP_focus");

    TS_ASSERT_EQUALS(outWS->getNumberHistograms(), 1);
    AnalysisDataService::Instance().remove("SNAP_focus");
  }

  void test_SNAP_event_six_groups_dontPreserveEvents() {
    auto alg = AlgorithmFactory::Instance().create("DiffractionFocussing", 2);
    alg->initialize();
    alg->setPropertyValue("InputWorkspace", "SNAP_empty");
    alg->setPropertyValue("GroupingWorkspace", "SNAP_group_several");
    alg->setPropertyValue("OutputWorkspace", "SNAP_focus");
    alg->setPropertyValue("PreserveEvents", "0");
    alg->execute();
    MatrixWorkspace_sptr outWS = AnalysisDataService::Instance().retrieveWS<MatrixWorkspace>("SNAP_focus");

    TS_ASSERT_EQUALS(outWS->getNumberHistograms(), 6);
    AnalysisDataService::Instance().remove("SNAP_focus");
  }
};
