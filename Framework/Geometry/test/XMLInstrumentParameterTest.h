// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidGeometry/Instrument/XMLInstrumentParameter.h"
#include "MantidKernel/Exception.h"
#include "MantidKernel/System.h"
#include "MantidKernel/TimeROI.h"
#include "MantidKernel/TimeSeriesProperty.h"
#include "MantidKernel/Timer.h"
#include <cxxtest/TestSuite.h>
#include <memory>

using namespace Mantid;
using namespace Mantid::Geometry;
using namespace Mantid::Kernel;

class XMLInstrumentParameterTest : public CxxTest::TestSuite {
private:
  using XMLInstrumentParameter_sptr = std::shared_ptr<XMLInstrumentParameter>;

  /**
  Construction logic for the XMLInstrumentParameter type isn't great, so this
  method acts a helper to keep the test methods cleaner.
  */
  XMLInstrumentParameter_sptr make_logfile_object(const std::string &filterBy) {
    const std::string logfileID = "1";
    const std::string value;
    const std::shared_ptr<Kernel::Interpolation> &interpolation = std::make_shared<Interpolation>();
    const std::string formula;
    const std::string formulaUnit;
    const std::string resultUnit;
    const std::string paramName;
    const std::string type;
    const std::string tie;
    const std::vector<std::string> constraint;
    std::string penaltyFactor;
    const std::string fitFunc;
    const std::string eq;
    const Geometry::IComponent *comp = nullptr;
    double angleConvertConst = 0.0;

    return std::shared_ptr<XMLInstrumentParameter>(new XMLInstrumentParameter(
        logfileID, value, interpolation, formula, formulaUnit, resultUnit, paramName, type, tie, constraint,
        penaltyFactor, fitFunc, filterBy, eq, comp, angleConvertConst, "", "true"));
  }

public:
  void test_throws_with_unknown_flag() {
    Kernel::TimeROI *roi = nullptr;
    TimeSeriesProperty<double> series("doubleProperty");
    series.addValue("2000-11-30T01:01:01", 1);

    const std::string made_up_flag = "mode"; // We do not support mode statistics filtering.
    XMLInstrumentParameter_sptr logFile = make_logfile_object(made_up_flag);

    TSM_ASSERT_THROWS("Unknown flag should cause failure", logFile->createParamValue(&series, roi),
                      const Kernel::Exception::InstrumentDefinitionError &)
  }

  void test_filter_by_first_value() {
    Kernel::TimeROI *roi = nullptr;
    TimeSeriesProperty<double> series("doubleProperty");

    const double expectedFilteredValue = 1;
    series.addValue("2000-11-30T01:01:01", expectedFilteredValue);
    series.addValue("2000-11-30T01:01:02", 2);

    XMLInstrumentParameter_sptr logFile = make_logfile_object("first_value");
    const double actualFilteredValue = logFile->createParamValue(&series, roi);
    TSM_ASSERT_EQUALS("Filtering by First Value is not performed correctly", expectedFilteredValue,
                      actualFilteredValue);
  }

  void test_filter_by_last_value() {
    Kernel::TimeROI *roi = new Kernel::TimeROI;
    roi->addROI("2000-11-30T01:01:01", "2000-11-30T01:01:02");
    TimeSeriesProperty<double> series("doubleProperty");

    const double expectedFilteredValue = 1;
    series.addValue("2000-11-30T01:01:01", 0);
    series.addValue("2000-11-30T01:01:02", expectedFilteredValue);

    XMLInstrumentParameter_sptr logFile = make_logfile_object("last_value");
    const double actualFilteredValue = logFile->createParamValue(&series, roi);
    TSM_ASSERT_EQUALS("Filtering by Last Value is not performed correctly", expectedFilteredValue, actualFilteredValue);
  }

  void test_filter_by_maximum_value() {
    Kernel::TimeROI *roi = new Kernel::TimeROI;
    roi->addROI("2000-11-30T01:01:01", "2000-11-30T01:01:03");
    TimeSeriesProperty<double> series("doubleProperty");

    const double expectedFilteredValue = 1;
    series.addValue("2000-11-30T01:01:01", 0.1);
    series.addValue("2000-11-30T01:01:02",
                    expectedFilteredValue); // maximum. 1 > 0.9 > 0.1
    series.addValue("2000-11-30T01:01:03", 0.9);

    XMLInstrumentParameter_sptr logFile = make_logfile_object("maximum");
    const double actualFilteredValue = logFile->createParamValue(&series, roi);
    TSM_ASSERT_EQUALS("Filtering by Maximum is not performed correctly", expectedFilteredValue, actualFilteredValue);
  }

  void test_filter_by_minimum_value() {
    Kernel::TimeROI *roi = new Kernel::TimeROI;
    roi->addROI("2000-11-30T01:01:01", "2000-11-30T01:01:03");
    TimeSeriesProperty<double> series("doubleProperty");

    const double expectedFilteredValue = 1;
    series.addValue("2000-11-30T01:01:01", 3);
    series.addValue("2000-11-30T01:01:02",
                    expectedFilteredValue); // minimum. 1 < 3 < 4
    series.addValue("2000-11-30T01:01:03", 4);

    XMLInstrumentParameter_sptr logFile = make_logfile_object("minimum");
    const double actualFilteredValue = logFile->createParamValue(&series, roi);
    TSM_ASSERT_EQUALS("Filtering by Minimum is not performed correctly", expectedFilteredValue, actualFilteredValue);
  }

  void test_filter_by_mean_value() {
    Kernel::TimeROI *roi = new Kernel::TimeROI;
    roi->addROI("2000-11-30T01:01:01", "2000-11-30T01:01:03");
    TimeSeriesProperty<double> series("doubleProperty");

    const double expectedFilteredValue = 1;
    series.addValue("2000-11-30T01:01:01", 0);
    series.addValue("2000-11-30T01:01:02",
                    expectedFilteredValue); // time series mean = value at T =
                                            // (T1 + T2 + T3) / 3
    series.addValue("2000-11-30T01:01:03", 2);

    XMLInstrumentParameter_sptr logFile = make_logfile_object("mean");
    const double actualFilteredValue = logFile->createParamValue(&series, roi);
    TSM_ASSERT_EQUALS("Filtering by Mean is not performed correctly", expectedFilteredValue, actualFilteredValue);
  }

  void test_filter_by_median_value() {
    Kernel::TimeROI *roi = new Kernel::TimeROI;
    roi->addROI("2000-11-30T01:01:01", "2000-11-30T01:02:00");
    TimeSeriesProperty<double> series("doubleProperty");

    const double expectedFilteredValue = 2;
    series.addValue("2000-11-30T01:01:01", 0);
    series.addValue("2000-11-30T01:01:02", 1);
    series.addValue("2000-11-30T01:01:03",
                    expectedFilteredValue); // Median time.
    series.addValue("2000-11-30T01:01:04", 4);
    series.addValue("2000-11-30T01:02:00", 5);

    XMLInstrumentParameter_sptr logFile = make_logfile_object("median");
    const double actualFilteredValue = logFile->createParamValue(&series, roi);
    TSM_ASSERT_EQUALS("Filtering by Median is not performed correctly", expectedFilteredValue, actualFilteredValue);
  }

  // This functionality will soon be legacy, since filtering by nth-position is
  // not a good idea.
  void test_filter_by_nth_position() {
    Kernel::TimeROI *roi = nullptr;
    TimeSeriesProperty<double> series("doubleProperty");

    const double expectedFilteredValue = 1;
    series.addValue("2000-11-30T01:01:01", 0);
    series.addValue("2000-11-30T01:01:02", expectedFilteredValue);
    series.addValue("2000-11-30T01:01:03", 2);
    series.addValue("2000-11-30T01:01:04", 3);

    XMLInstrumentParameter_sptr logFile = make_logfile_object("position 2");
    const double actualFilteredValue = logFile->createParamValue(&series, roi);
    TSM_ASSERT_EQUALS("Filtering by Nth position is not performed correctly", expectedFilteredValue,
                      actualFilteredValue);
  }

  void test_help_string() {
    Kernel::TimeROI *roi = nullptr;
    std::vector<std::string> constr;
    std::string penaltyFactor;
    XMLInstrumentParameter testPar("logfileID", "value", std::make_shared<Interpolation>(), "formula", "sourceFU",
                                   "resultFU", "testPar", "aType", "noTie", constr, penaltyFactor, "aFitFunc",
                                   "FilterBy", "eqTo", nullptr, 0.0, "test string.     Long test string.", "true");

    TS_ASSERT_EQUALS(testPar.m_description, "test string. Long test string.");
  }

  void test_parameter_not_visible() {
    Kernel::TimeROI *roi = nullptr;
    std::vector<std::string> constr;
    std::string penaltyFactor;
    XMLInstrumentParameter testPar("logfileID", "value", std::make_shared<Interpolation>(), "formula", "sourceFU",
                                   "resultFU", "testPar", "aType", "noTie", constr, penaltyFactor, "aFitFunc",
                                   "FilterBy", "eqTo", nullptr, 0.0, "test string.     Long test string.", "false");

    TS_ASSERT_EQUALS(testPar.m_visible, "false");
  }
};
