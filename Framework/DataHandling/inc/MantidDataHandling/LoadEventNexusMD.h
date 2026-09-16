// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidAPI/BoxControllerSettingsAlgorithm.h"
#include "MantidAPI/IMDEventWorkspace_fwd.h"
#include "MantidAPI/MatrixWorkspace_fwd.h"
#include "MantidDataHandling/DllConfig.h"
#include "MantidDataHandling/LoadEventNexusMD/QConversionTable.h"
#include "MantidGeometry/MDGeometry/MDTypes.h"
#include "MantidNexus/NexusDescriptor.h"

#include <H5Cpp.h>
#include <atomic>
#include <string>
#include <utility>
#include <vector>

namespace Mantid::DataHandling::LoadEventNexusMD {

/**
 * Load single crystal diffraction event data straight from a NeXus file into an MDEventWorkspace in Q_sample or
 * Q_lab, without going through an EventWorkspace. Equivalent to LoadEventNexus followed by ConvertToMD, but the
 * per-pixel event lists are never built.
 */
class MANTID_DATAHANDLING_DLL LoadEventNexusMD : public API::BoxControllerSettingsAlgorithm {
public:
  const std::string name() const override;
  int version() const override;
  const std::string category() const override;
  const std::string summary() const override;
  const std::vector<std::string> seeAlso() const override;

private:
  void init() override;
  std::map<std::string, std::string> validateInputs() override;
  void exec() override;

  /// Collect the names of the NXevent_data entries in the file
  std::vector<std::string> determineBanksToLoad(const Mantid::Nexus::NexusDescriptor &descriptor) const;
  /// Additive time-of-flight correction in microseconds, from the instrument's T0 parameter
  double getT0(const API::MatrixWorkspace_sptr &wksp) const;
  /// Smallest and largest corrected time-of-flight in the file, in microseconds. Reads event_time_offset only.
  std::pair<double, double> findTimeOfFlightRange(H5::H5File &h5file, const std::vector<std::string> &bankEntryNames,
                                                  const double t0) const;
  /// The box enclosing every event, from the per-detector directions and the time-of-flight range
  void calculateExtents(H5::H5File &h5file, const std::vector<std::string> &bankEntryNames,
                        const QConversionTable &table, const double t0, std::vector<double> &minVals,
                        std::vector<double> &maxVals) const;
  /// Read every bank in chunks, convert to Q and add to the box tree, splitting between chunks
  template <typename MDEW>
  void loadEvents(const std::shared_ptr<MDEW> &ws, H5::H5File &h5file, const std::vector<std::string> &bankEntryNames,
                  const QConversionTable &table, const double t0, const Mantid::coord_t *extentsMin,
                  const Mantid::coord_t *extentsMax);

  /// events whose event_id matches no detector in the instrument definition, as LoadEventNexus also discards
  std::atomic<size_t> m_eventsWithNoDetector{0};
};

// these properties are public to simplify testing and calling from other code
namespace PropertyNames {
const std::string FILENAME("Filename");
const std::string OUTPUT_WKSP("OutputWorkspace");
const std::string Q_FRAME("QFrame");
const std::string EVENT_TYPE("OutputEventType");
const std::string MIN_VALUES("MinValues");
const std::string MAX_VALUES("MaxValues");
const std::string READ_SIZE_FROM_DISK("ReadSizeFromDisk");
const std::string EVENTS_PER_THREAD("EventsPerThread");
} // namespace PropertyNames

namespace QFrames {
const std::string SAMPLE("Q (sample frame)");
const std::string LAB("Q (lab frame)");
} // namespace QFrames

namespace EventTypes {
const std::string FULL("MDEvent");
const std::string LEAN("MDLeanEvent");
} // namespace EventTypes

} // namespace Mantid::DataHandling::LoadEventNexusMD
