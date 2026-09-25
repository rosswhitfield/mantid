// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "MantidDataHandling/AlignAndFocusPowderSlim/ProcessFileOrderTask.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/ProcessEventsTask.h"
#include "MantidKernel/Logger.h"
#include "MantidNexus/H5Util.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <deque>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

namespace {
auto g_log = Kernel::Logger("ProcessFileOrderTask");

/// Events read and held per wave, in bytes
constexpr uint64_t WAVE_BYTES = 256 * 1024 * 1024;
/// Wave buffers: one wave is histogrammed while the reads of the next NUM_BUFFERS - 1 are queued, so the reading
/// threads always have work at the end of a wave instead of idling until the next one is started
constexpr size_t NUM_BUFFERS = 4;
/// Largest single read, which is also the size of each reading thread's staging buffer
constexpr uint64_t MAX_SPAN = 8 * 1024 * 1024;
constexpr uint64_t EVENT_BYTES = sizeof(uint32_t) + sizeof(float);

double secondsSince(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

std::string columnPath(const std::string &bank, const std::string &column) { return "/entry/" + bank + "/" + column; }
} // namespace

bool ProcessFileOrderTask::canProcess(const DirectEventReader &reader, const std::vector<std::string> &bankEntryNames) {
  return std::all_of(bankEntryNames.cbegin(), bankEntryNames.cend(), [&reader](const std::string &bank) {
    return bank.empty() || (reader.column(columnPath(bank, NxsFieldNames::DETID)) != nullptr &&
                            reader.column(columnPath(bank, NxsFieldNames::TIME_OF_FLIGHT)) != nullptr);
  });
}

ProcessFileOrderTask::ProcessFileOrderTask(std::vector<std::string> &bankEntryNames, H5::H5File &h5file,
                                           std::shared_ptr<NexusLoader> loader,
                                           std::shared_ptr<const DirectEventReader> reader,
                                           SpectraProcessingData &processingData,
                                           const BankCalibrationFactory &calibFactory, const size_t events_per_chunk,
                                           const size_t grainsize_event)
    : ProcessBankTaskBase(bankEntryNames, std::move(loader), calibFactory), m_reader(std::move(reader)),
      m_processingData(processingData), m_grainsize_event(grainsize_event) {
  auto entry = h5file.openGroup("entry");
  for (size_t bank_index = 0; bank_index < bankEntryNames.size(); ++bank_index) {
    const auto &bankName = this->bankName(bank_index);
    if (bankName.empty())
      continue;
    auto event_group = entry.openGroup(bankName);
    auto tof_SDS = event_group.openDataSet(NxsFieldNames::TIME_OF_FLIGHT);
    const auto total_events = static_cast<uint64_t>(tof_SDS.getSpace().getSelectNpoints());
    if (total_events == 0)
      continue;
    std::string tof_unit;
    Nexus::H5Util::readStringAttribute(tof_SDS, "units", tof_unit);
    auto calibrations = this->getCalibrations(tof_unit, bank_index);
    if (calibrations.empty()) {
      g_log.debug() << "skipping " << bankName << " because calibration is empty\n";
      continue;
    }

    const auto detidPath = columnPath(bankName, NxsFieldNames::DETID);
    const auto tofPath = columnPath(bankName, NxsFieldNames::TIME_OF_FLIGHT);
    const auto *detidColumn = m_reader->column(detidPath);
    const auto *tofColumn = m_reader->column(tofPath);
    if (detidColumn == nullptr || tofColumn == nullptr)
      throw std::logic_error(bankName + " cannot be read directly; check canProcess() first");
    const size_t bank = m_banks.size();
    m_banks.push_back({detidPath, tofPath, std::move(calibrations)});

    // cut the bank's events into batches of at most events_per_chunk, the same way ProcessBankTask does
    auto eventRanges = this->getEventIndexRanges(event_group, total_events);
    while (!eventRanges.empty()) {
      Batch batch;
      batch.bank = bank;
      while (!eventRanges.empty() && batch.numEvents < events_per_chunk) {
        const auto eventRange = eventRanges.top();
        eventRanges.pop();
        const size_t range_size = eventRange.second - eventRange.first;
        const size_t remaining = events_per_chunk - batch.numEvents;
        const size_t take = std::min(range_size, remaining);
        batch.offsets.push_back(eventRange.first);
        batch.slabsizes.push_back(take);
        batch.numEvents += take;
        if (take < range_size) {
          eventRanges.emplace(eventRange.first + take, eventRange.second);
          break;
        }
      }
      if (batch.numEvents == 0)
        continue;
      // roughly where the batch starts in the file, to order the batches of all banks by position; the chunks are
      // located later, as their waves come up
      const auto first = batch.offsets.front();
      batch.filePosition = std::min(m_reader->positionHint(detidPath, first), m_reader->positionHint(tofPath, first));
      m_batches.push_back(std::move(batch));
    }
  }

  std::stable_sort(m_batches.begin(), m_batches.end(),
                   [](const Batch &left, const Batch &right) { return left.filePosition < right.filePosition; });
  size_t first = 0;
  uint64_t bytes = 0;
  for (size_t i = 0; i < m_batches.size(); ++i) {
    bytes += m_batches[i].numEvents * EVENT_BYTES;
    if (bytes >= WAVE_BYTES || i + 1 == m_batches.size()) {
      m_waves.emplace_back(first, i + 1);
      m_waveCapacity = std::max<size_t>(m_waveCapacity, bytes / EVENT_BYTES);
      first = i + 1;
      bytes = 0;
    }
  }
  g_log.debug() << m_batches.size() << " batches from " << m_banks.size() << " banks in " << m_waves.size()
                << " waves\n";
}

std::vector<ParallelFileReader::Task> ProcessFileOrderTask::prepareWave(const Wave &wave, WaveBuffer &buffer) {
  // locate the wave's chunks first, reading the index leaves that hold them in one round
  const auto locateStart = std::chrono::steady_clock::now();
  std::vector<ElementRange> ranges;
  for (size_t index = wave.first; index < wave.second; ++index) {
    const auto &batch = m_batches[index];
    const auto &bank = m_banks[batch.bank];
    for (size_t slab = 0; slab < batch.offsets.size(); ++slab) {
      ranges.push_back({bank.detidPath, batch.offsets[slab], batch.slabsizes[slab]});
      ranges.push_back({bank.tofPath, batch.offsets[slab], batch.slabsizes[slab]});
    }
  }
  m_reader->locateChunks(ranges);
  m_timing.findingChunks += secondsSince(locateStart);

  std::vector<ByteRun> runs;
  size_t used = 0;
  for (size_t index = wave.first; index < wave.second; ++index) {
    auto &batch = m_batches[index];
    batch.detid = std::span<uint32_t>(buffer.detid.get() + used, batch.numEvents);
    batch.tof = std::span<float>(buffer.tof.get() + used, batch.numEvents);
    used += batch.numEvents;
    const auto &bank = m_banks[batch.bank];
    const auto *detidColumn = m_reader->column(bank.detidPath);
    const auto *tofColumn = m_reader->column(bank.tofPath);
    char *detidDest = reinterpret_cast<char *>(batch.detid.data());
    char *tofDest = reinterpret_cast<char *>(batch.tof.data());
    for (size_t slab = 0; slab < batch.offsets.size(); ++slab) {
      for (const auto &run : planRuns(*detidColumn, batch.offsets[slab], batch.slabsizes[slab], detidDest))
        runs.push_back(run);
      for (const auto &run : planRuns(*tofColumn, batch.offsets[slab], batch.slabsizes[slab], tofDest))
        runs.push_back(run);
      detidDest += batch.slabsizes[slab] * sizeof(uint32_t);
      tofDest += batch.slabsizes[slab] * sizeof(float);
    }
  }

  // reading through gaps of up to one block costs nothing, since reads are widened to whole blocks anyway
  const auto block = m_reader->blockSize();
  auto spans = planSpans(std::move(runs), block, MAX_SPAN, block, m_reader->fileReader().fileSize());

  std::vector<ParallelFileReader::Task> tasks;
  tasks.reserve(spans.size());
  for (auto &span : spans) {
    tasks.emplace_back([span = std::move(span)](std::istream &stream) {
      thread_local std::vector<char> staging;
      if (staging.size() < span.size)
        staging.resize(span.size);
      readAt(stream, span.fileOffset, span.size, staging.data());
      for (const auto &piece : span.pieces)
        std::memcpy(piece.dest, staging.data() + (piece.fileOffset - span.fileOffset), piece.size);
    });
  }
  return tasks;
}

void ProcessFileOrderTask::histogramWave(const Wave &wave) {
  tbb::parallel_for(
      tbb::blocked_range<size_t>(wave.first, wave.second, 1), [&](const tbb::blocked_range<size_t> &batch_range) {
        for (size_t index = batch_range.begin(); index < batch_range.end(); ++index) {
          auto &batch = m_batches[index];
          const auto &calibrations = m_banks[batch.bank].calibrations;
          tbb::parallel_for(tbb::blocked_range<size_t>(0, m_processingData.counts.size()),
                            [&](const tbb::blocked_range<size_t> &output_range) {
                              for (size_t output_index = output_range.begin(); output_index < output_range.end();
                                   ++output_index) {
                                ProcessEventsTask task(&batch.detid, &batch.tof, &calibrations.at(output_index),
                                                       m_processingData.binedges[output_index]);
                                const tbb::blocked_range<size_t> range_info(0, batch.numEvents, m_grainsize_event);
                                tbb::parallel_reduce(range_info, task);
                                auto &counts = m_processingData.counts[output_index];
                                for (size_t i = 0; i < counts.size(); ++i)
                                  counts[i].fetch_add(task.y_temp[i], std::memory_order_relaxed);
                              }
                            });
        }
      });
}

void ProcessFileOrderTask::run(API::Progress &progress) {
  if (m_waves.empty())
    return;
  // allocated without initialising: every element is written by a read before it is used
  const size_t numBuffers = std::min(NUM_BUFFERS, m_waves.size());
  std::vector<WaveBuffer> buffers(numBuffers);
  for (auto &buffer : buffers) {
    buffer.detid.reset(new uint32_t[m_waveCapacity]);
    buffer.tof.reset(new float[m_waveCapacity]);
  }
  auto &fileReader = m_reader->fileReader();
  std::deque<std::future<void>> pending; // the reads of waves [index, next), in order
  size_t next = 0;
  try {
    for (size_t index = 0; index < m_waves.size(); ++index) {
      // Queue waves while a buffer is free. Wave w uses buffer w % numBuffers, last used by wave w - numBuffers, which
      // has been histogrammed once it is below index.
      auto start = std::chrono::steady_clock::now();
      while (next < m_waves.size() && next < index + numBuffers) {
        pending.push_back(fileReader.start(prepareWave(m_waves[next], buffers[next % numBuffers])));
        ++next;
      }
      m_timing.preparingReads += secondsSince(start);

      start = std::chrono::steady_clock::now();
      auto current = std::move(pending.front());
      pending.pop_front();
      current.get();
      m_timing.waitingForReads += secondsSince(start);

      start = std::chrono::steady_clock::now();
      histogramWave(m_waves[index]);
      m_timing.histogramming += secondsSince(start);
      progress.report();
    }
  } catch (...) {
    // queued reads write into the buffers; let them finish before the buffers go away
    for (auto &reads : pending) {
      if (reads.valid())
        reads.wait();
    }
    throw;
  }
}

} // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
