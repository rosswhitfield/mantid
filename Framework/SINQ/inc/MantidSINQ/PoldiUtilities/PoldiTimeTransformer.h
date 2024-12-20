// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2014 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidSINQ/DllConfig.h"
#include "MantidSINQ/PoldiUtilities/PoldiConversions.h"
#include "MantidSINQ/PoldiUtilities/PoldiInstrumentAdapter.h"

namespace Mantid {
namespace Poldi {

/** PoldiTimeTransformer

    This helper class transforms peaks from "d" to "time", using
    factors derived from POLDI detector configuration.

      @author Michael Wedel, Paul Scherrer Institut - SINQ
      @date 22/05/2014
  */

struct DetectorElementCharacteristics {
  DetectorElementCharacteristics()
      : distance(0.0), totalDistance(0.0), twoTheta(0.0), sinTheta(0.0), cosTheta(1.0), tof1A(0.0) {}

  DetectorElementCharacteristics(int element, const PoldiAbstractDetector_sptr &detector,
                                 const PoldiAbstractChopper_sptr &chopper) {
    distance = detector->distanceFromSample(element);
    totalDistance = detector->distanceFromSample(element) + chopper->distanceFromSample();
    twoTheta = detector->twoTheta(element);
    sinTheta = sin(twoTheta / 2.0);
    cosTheta = cos(twoTheta / 2.0);
    tof1A = Conversions::dtoTOF(1.0, totalDistance, sinTheta);
  }

  double distance;
  double totalDistance;
  double twoTheta;
  double sinTheta;
  double cosTheta;
  double tof1A;
};

class DetectorElementData {
public:
  DetectorElementData(int element, const DetectorElementCharacteristics &center,
                      const PoldiAbstractDetector_sptr &detector, const PoldiAbstractChopper_sptr &chopper) {
    DetectorElementCharacteristics current(element, detector, chopper);

    m_intensityFactor = pow(center.distance / current.distance, 2.0) * current.sinTheta / center.sinTheta;
    m_lambdaFactor = 2.0 * current.sinTheta / center.tof1A;
    m_timeFactor = current.sinTheta / center.sinTheta * current.totalDistance / center.totalDistance;
    m_widthFactor = current.cosTheta - center.cosTheta;
    m_tofFactor = center.tof1A / current.tof1A;
  }

  double intensityFactor() const { return m_intensityFactor; }
  double lambdaFactor() const { return m_lambdaFactor; }
  double timeFactor() const { return m_timeFactor; }
  double widthFactor() const { return m_widthFactor; }
  double tofFactor() const { return m_tofFactor; }

protected:
  double m_intensityFactor;
  double m_lambdaFactor;
  double m_timeFactor;
  double m_widthFactor;
  double m_tofFactor;
};

using DetectorElementData_const_sptr = std::shared_ptr<const DetectorElementData>;

class MANTID_SINQ_DLL PoldiTimeTransformer {
public:
  PoldiTimeTransformer();
  PoldiTimeTransformer(const PoldiInstrumentAdapter_sptr &poldiInstrument);
  virtual ~PoldiTimeTransformer() = default;

  void initializeFromPoldiInstrument(const PoldiInstrumentAdapter_sptr &poldiInstrument);

  size_t detectorElementCount() const;

  double dToTOF(double d) const;
  double detectorElementIntensity(double centreD, size_t detectorIndex) const;

  double calculatedTotalIntensity(double centreD) const;

protected:
  std::vector<DetectorElementData_const_sptr> getDetectorElementData(const PoldiAbstractDetector_sptr &detector,
                                                                     const PoldiAbstractChopper_sptr &chopper);
  DetectorElementCharacteristics getDetectorCenterCharacteristics(const PoldiAbstractDetector_sptr &detector,
                                                                  const PoldiAbstractChopper_sptr &chopper);

  DetectorElementCharacteristics m_detectorCenter;
  std::vector<DetectorElementData_const_sptr> m_detectorElementData;
  double m_detectorEfficiency;
  size_t m_chopperSlits;

  PoldiSourceSpectrum_const_sptr m_spectrum;
};

using PoldiTimeTransformer_sptr = std::shared_ptr<PoldiTimeTransformer>;

} // namespace Poldi
} // namespace Mantid
