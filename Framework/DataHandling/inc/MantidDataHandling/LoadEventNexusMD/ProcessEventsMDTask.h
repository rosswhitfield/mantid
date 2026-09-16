// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#pragma once

#include "MantidDataHandling/LoadEventNexusMD/QConversionTable.h"
#include "MantidDataObjects/MDEventInserter.h"
#include "MantidGeometry/IDTypes.h"
#include "MantidGeometry/MDGeometry/MDTypes.h"

#include <atomic>
#include <memory>
#include <tbb/blocked_range.h>
#include <vector>

namespace Mantid::DataHandling::LoadEventNexusMD {

/**
 * tbb::parallel_for body that turns a chunk of (detector id, time-of-flight) pairs from one bank into MDEvents and
 * adds them to the output workspace.
 *
 * Events are added one at a time rather than through MDEventWorkspace::addEvents. MDGridBox does not override
 * addEvents, so the bulk call takes the mutex of the *root* box and serialises every thread on a single lock, whereas
 * MDGridBox::addEvent recurses down to a leaf MDBox whose mutex is per-leaf and therefore scales.
 *
 * The caller must not split the box tree while any of these are running.
 */
template <typename MDEW> class ProcessEventsMDTask {
public:
  /**
   * @param tofScale Factor converting this bank's time-of-flight into microseconds.
   * @param tofOffset Added after scaling, the instrument's T0 parameter.
   */
  ProcessEventsMDTask(const std::vector<uint32_t> *detids, const std::vector<float> *tofs,
                      const QConversionTable *table, std::shared_ptr<MDEW> ws, const Mantid::coord_t *extentsMin,
                      const Mantid::coord_t *extentsMax, const double tofScale, const double tofOffset,
                      std::atomic<size_t> *eventsAdded, std::atomic<size_t> *eventsWithNoDetector)
      : m_detids(detids), m_tofs(tofs), m_table(table), m_ws(std::move(ws)), m_extentsMin(extentsMin),
        m_extentsMax(extentsMax), m_tofScale(tofScale), m_tofOffset(tofOffset), m_eventsAdded(eventsAdded),
        m_eventsWithNoDetector(eventsWithNoDetector) {}

  void operator()(const tbb::blocked_range<size_t> &range) const {
    auto ws = m_ws; // MDEventInserter takes a non-const reference
    DataObjects::MDEventInserter<std::shared_ptr<MDEW>> inserter(ws);

    size_t added{0};
    size_t noDetector{0};
    for (size_t i = range.begin(); i < range.end(); ++i) {
      const auto detid = static_cast<detid_t>((*m_detids)[i]);
      const auto &coeff = m_table->value(detid);
      if (coeff.kfactor <= 0.) {
        ++noDetector; // a monitor, or an event_id the instrument definition does not describe
        continue;
      }

      // |k| = 2*pi/lambda. A zero or negative time-of-flight gives a non-finite coordinate, which the extents check
      // below rejects.
      const double tof = static_cast<double>((*m_tofs)[i]) * m_tofScale + m_tofOffset;
      const double k = coeff.kfactor / tof;
      Mantid::coord_t center[3] = {static_cast<Mantid::coord_t>(coeff.qx * k),
                                   static_cast<Mantid::coord_t>(coeff.qy * k),
                                   static_cast<Mantid::coord_t>(coeff.qz * k)};

      if (center[0] < m_extentsMin[0] || !(center[0] < m_extentsMax[0]))
        continue;
      if (center[1] < m_extentsMin[1] || !(center[1] < m_extentsMax[1]))
        continue;
      if (center[2] < m_extentsMin[2] || !(center[2] < m_extentsMax[2]))
        continue;

      inserter.insertMDEvent(1.f, 1.f, 0, 0, static_cast<int32_t>(detid), center);
      ++added;
    }

    if (added > 0)
      m_eventsAdded->fetch_add(added, std::memory_order_relaxed);
    if (noDetector > 0)
      m_eventsWithNoDetector->fetch_add(noDetector, std::memory_order_relaxed);
  }

private:
  const std::vector<uint32_t> *m_detids;
  const std::vector<float> *m_tofs;
  const QConversionTable *m_table;
  std::shared_ptr<MDEW> m_ws;
  const Mantid::coord_t *m_extentsMin;
  const Mantid::coord_t *m_extentsMax;
  double m_tofScale;
  double m_tofOffset;
  std::atomic<size_t> *m_eventsAdded;
  std::atomic<size_t> *m_eventsWithNoDetector;
};

} // namespace Mantid::DataHandling::LoadEventNexusMD
