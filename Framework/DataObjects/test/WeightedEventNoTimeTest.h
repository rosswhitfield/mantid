// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidDataObjects/Events.h"
#include "MantidKernel/Timer.h"
#include <cmath>
#include <cxxtest/TestSuite.h>

using namespace Mantid;
using namespace Mantid::Kernel;
using namespace Mantid::DataObjects;
using Mantid::Types::Event::TofEvent;

using std::runtime_error;
using std::size_t;
using std::vector;

//==========================================================================================
class WeightedEventNoTimeTest : public CxxTest::TestSuite {
public:
  void testConstructors() {
    TofEvent e(123, 456);
    WeightedEvent we;
    WeightedEventNoTime wen, wen2;

    // Empty
    wen = WeightedEventNoTime();
    TS_ASSERT_EQUALS(wen.tof(), 0);
    TS_ASSERT_EQUALS(wen.pulseTime(), 0);
    TS_ASSERT_EQUALS(wen.weight(), 1.0);
    TS_ASSERT_EQUALS(wen.error(), 1.0);

    // From WeightedEvent
    we = WeightedEvent(456, 789, 2.5, 1.5 * 1.5);
    wen = WeightedEventNoTime(we);
    TS_ASSERT_EQUALS(wen.tof(), 456);
    TS_ASSERT_EQUALS(wen.pulseTime(), 0); // Lost the time!
    TS_ASSERT_EQUALS(wen.weight(), 2.5);
    TS_ASSERT_EQUALS(wen.error(), 1.5);

    // Default one weight from TofEvent
    wen = WeightedEventNoTime(e);
    TS_ASSERT_EQUALS(wen.tof(), 123);
    TS_ASSERT_EQUALS(wen.pulseTime(), 0);
    TS_ASSERT_EQUALS(wen.weight(), 1.0);
    TS_ASSERT_EQUALS(wen.error(), 1.0);

    // TofEvent + weights
    wen = WeightedEventNoTime(e, 3.5, 0.5 * 0.5);
    TS_ASSERT_EQUALS(wen.tof(), 123);
    TS_ASSERT_EQUALS(wen.pulseTime(), 0);
    TS_ASSERT_EQUALS(wen.weight(), 3.5);
    TS_ASSERT_EQUALS(wen.error(), 0.5);

    // Full constructor
    wen = WeightedEventNoTime(456, 2.5, 1.5 * 1.5);
    TS_ASSERT_EQUALS(wen.tof(), 456);
    TS_ASSERT_EQUALS(wen.pulseTime(), 0); // Never had time
    TS_ASSERT_EQUALS(wen.weight(), 2.5);
    TS_ASSERT_EQUALS(wen.error(), 1.5);
  }

  void testAssignAndCopy_AndEquality() {
    WeightedEventNoTime wen, wen2;

    // Copy constructor
    wen = WeightedEventNoTime();
    wen2 = WeightedEventNoTime(456, 2.5, 1.5 * 1.5);
    TS_ASSERT(!(wen == wen2));

    wen = wen2;
    TS_ASSERT_EQUALS(wen.tof(), 456);
    TS_ASSERT_EQUALS(wen.pulseTime(), 0); // Never had time
    TS_ASSERT_EQUALS(wen.weight(), 2.5);
    TS_ASSERT_EQUALS(wen.error(), 1.5);

    TS_ASSERT(wen == wen2);
  }
};
