// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidSINQ/PoldiUtilities/PoldiDeadWireDecorator.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"

#include <algorithm>
#include <utility>

namespace Mantid::Poldi {

using namespace Geometry;

PoldiDeadWireDecorator::PoldiDeadWireDecorator(std::set<int> deadWires,
                                               const std::shared_ptr<Poldi::PoldiAbstractDetector> &detector)
    : PoldiDetectorDecorator(detector), m_deadWireSet(std::move(deadWires)), m_goodElements() {
  setDecoratedDetector(detector);
}

PoldiDeadWireDecorator::PoldiDeadWireDecorator(const Geometry::DetectorInfo &poldiDetectorInfo,
                                               const std::shared_ptr<PoldiAbstractDetector> &detector)
    : PoldiDetectorDecorator(detector), m_deadWireSet(), m_goodElements() {
  setDecoratedDetector(detector);

  std::vector<detid_t> allDetectorIds = poldiDetectorInfo.detectorIDs();
  std::vector<detid_t> deadDetectorIds(allDetectorIds.size());

  auto endIterator = std::remove_copy_if(
      allDetectorIds.begin(), allDetectorIds.end(), deadDetectorIds.begin(),
      [&](const detid_t detID) -> bool { return !poldiDetectorInfo.isMasked(poldiDetectorInfo.indexOf(detID)); });
  deadDetectorIds.resize(std::distance(deadDetectorIds.begin(), endIterator));

  setDeadWires(std::set<int>(deadDetectorIds.begin(), deadDetectorIds.end()));
}

void PoldiDeadWireDecorator::setDeadWires(std::set<int> deadWires) {
  m_deadWireSet = std::move(deadWires);

  detectorSetHook();
}

const std::set<int> &PoldiDeadWireDecorator::deadWires() { return m_deadWireSet; }

size_t PoldiDeadWireDecorator::elementCount() { return m_goodElements.size(); }

const std::vector<int> &PoldiDeadWireDecorator::availableElements() { return m_goodElements; }

void PoldiDeadWireDecorator::detectorSetHook() {
  if (m_decoratedDetector) {
    m_goodElements = getGoodElements(m_decoratedDetector->availableElements());
  } else {
    throw(std::runtime_error("No decorated detector set!"));
  }
}

std::vector<int> PoldiDeadWireDecorator::getGoodElements(std::vector<int> rawElements) {
  if (!m_deadWireSet.empty()) {
    if (*m_deadWireSet.rbegin() > rawElements.back()) {
      throw std::runtime_error(std::string("Deadwires set contains illegal index."));
    }
    size_t newElementCount = rawElements.size() - m_deadWireSet.size();

    std::vector<int> goodElements(newElementCount);
    using namespace std::placeholders;
    std::remove_copy_if(rawElements.begin(), rawElements.end(), goodElements.begin(),
                        std::bind(&PoldiDeadWireDecorator::isDeadElement, this, _1));

    return goodElements;
  }

  return rawElements;
}

bool PoldiDeadWireDecorator::isDeadElement(int index) { return m_deadWireSet.count(index) != 0; }

} // namespace Mantid::Poldi
