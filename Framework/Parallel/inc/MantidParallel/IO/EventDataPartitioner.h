#ifndef MANTID_PARALLEL_EVENTDATAPARTITIONER_H_
#define MANTID_PARALLEL_EVENTDATAPARTITIONER_H_

#include "MantidParallel/DllConfig.h"
#include "MantidParallel/IO/Chunker.h"
#include "MantidParallel/IO/PulseTimeGenerator.h"
#include "MantidTypes/Core/DateAndTime.h"

namespace Mantid {
namespace Parallel {
namespace IO {

/** Partition the event_time_offset and event_id entries and combine them with
  pulse time information obtained from PulseTimeGenerator. Partitioning is to
  obtain a separate vector of events for each rank in an MPI run of Mantid,
  i.e., each event_id is assigned to a specific MPI rank. Currently a
  round-robin partitioning scheme is hard-coded.

  @author Simon Heybrock
  @date 2017

  Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
  National Laboratory & European Spallation Source

  This file is part of Mantid.

  Mantid is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  Mantid is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  File change history is stored at: <https://github.com/mantidproject/mantid>
  Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
namespace detail {
template <class TimeOffsetType> struct Event {
  int32_t index; // local spectrum index
  TimeOffsetType tof;
  Types::Core::DateAndTime pulseTime;
};
}

template <class TimeOffsetType> class AbstractEventDataPartitioner {
public:
  using Event = detail::Event<TimeOffsetType>;
  AbstractEventDataPartitioner(const int numWorkers) : m_numWorkers(numWorkers) {}
  virtual ~AbstractEventDataPartitioner() = default;

  /** Partition given data.
   *
   * @param rankData output vector data vectors for each MPI rank.
   * @param globalSpectrumIndex list of spectrum indices
   * @param eventTimeOffset list TOF values, same length as globalSpectrumIndex
   * @param range defines start and end of data for lookup in PulseTimeGenerator
   */
  virtual void partition(std::vector<std::vector<Event>> &rankData,
                         const int32_t *globalSpectrumIndex,
                         const TimeOffsetType *eventTimeOffset,
                         const Chunker::LoadRange &range) = 0;

protected:
  const int m_numWorkers;
};

template <class IndexType, class TimeZeroType, class TimeOffsetType>
class EventDataPartitioner
    : public AbstractEventDataPartitioner<TimeOffsetType> {
public:
  using Event = detail::Event<TimeOffsetType>;
  EventDataPartitioner(const int numWorkers,
                       PulseTimeGenerator<IndexType, TimeZeroType> &&gen)
      : AbstractEventDataPartitioner<TimeOffsetType>(numWorkers),
        m_pulseTimes(std::move(gen)) {}

  void partition(std::vector<std::vector<Event>> &rankData,
                 const int32_t *globalSpectrumIndex,
                 const TimeOffsetType *eventTimeOffset,
                 const Chunker::LoadRange &range) override;

private:
  PulseTimeGenerator<IndexType, TimeZeroType> m_pulseTimes;
};

template <class IndexType, class TimeZeroType, class TimeOffsetType>
void EventDataPartitioner<
    IndexType, TimeZeroType,
    TimeOffsetType>::partition(std::vector<std::vector<Event>> &rankData,
                               const int32_t *globalSpectrumIndex,
                               const TimeOffsetType *eventTimeOffset,
                               const Chunker::LoadRange &range) {
  for (auto &item : rankData)
    item.clear();
  const auto workers =
      AbstractEventDataPartitioner<TimeOffsetType>::m_numWorkers;
  rankData.resize(workers);

  m_pulseTimes.seek(range.eventOffset);
  for (size_t event = 0; event < range.eventCount; ++event) {
    // Currently this supports only a hard-coded round-robin partitioning.
    int rank = globalSpectrumIndex[event] % workers;
    auto index = globalSpectrumIndex[event] / workers;
    rankData[rank].emplace_back(detail::Event<TimeOffsetType>{
        index, eventTimeOffset[event], m_pulseTimes.next()});
  }
}

} // namespace IO
} // namespace Parallel
} // namespace Mantid

#endif /* MANTID_PARALLEL_EVENTDATAPARTITIONER_H_ */
