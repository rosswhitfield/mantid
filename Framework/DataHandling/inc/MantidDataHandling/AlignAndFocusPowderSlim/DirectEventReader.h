// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#pragma once

#include "MantidDataHandling/DllConfig.h"

#include <H5Cpp.h>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <functional>
#include <future>
#include <istream>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

/** A pool of threads that read one file, each through its own stream.
 *
 * Reads are blocking system calls, so they run on these threads rather than
 * as TBB tasks: a TBB thread waiting on the disk can hold up unrelated work it
 * has stolen, and the file system needs many requests in flight to deliver
 * its bandwidth, far more than there are cores to spare.
 */
class MANTID_DATAHANDLING_DLL ParallelFileReader {
public:
  using Task = std::function<void(std::istream &)>;

  ParallelFileReader(const std::string &filename, size_t numThreads);
  ~ParallelFileReader();
  ParallelFileReader(const ParallelFileReader &) = delete;
  ParallelFileReader &operator=(const ParallelFileReader &) = delete;

  /// Start the tasks, each on a pool thread with that thread's stream. The future is ready when all have finished and
  /// carries the first exception thrown. Urgent tasks go ahead of everything already queued.
  std::future<void> start(std::vector<Task> tasks, bool urgent = false);
  /// Run the tasks and wait for them, rethrowing the first exception thrown.
  void run(std::vector<Task> tasks, bool urgent = false);

  uint64_t fileSize() const { return m_fileSize; }
  const std::string &filename() const { return m_filename; }

private:
  void work();

  const std::string m_filename;
  uint64_t m_fileSize;
  std::mutex m_mutex;
  std::condition_variable m_wake;
  std::deque<std::function<void(std::istream &)>> m_queue;
  bool m_stop{false};
  std::vector<std::thread> m_threads;
};

/// Read exactly @p size bytes at @p offset into @p dest; throws on a short read.
MANTID_DATAHANDLING_DLL void readAt(std::istream &stream, uint64_t offset, uint64_t size, char *dest);

/// Where each chunk of an unfiltered, chunked, one-dimensional dataset is stored in the file.
struct MANTID_DATAHANDLING_DLL ChunkedColumn {
  uint64_t numElements{0};
  uint64_t chunkElements{0};
  uint32_t elementSize{0};
  /// File offset of each chunk, in element order. Filled in as the chunks are located: see
  /// DirectEventReader::locateChunks.
  std::vector<uint64_t> chunkOffsets;
};

/// Elements [first, first + count) of the column at @p path.
struct ElementRange {
  std::string path;
  uint64_t first;
  uint64_t count;
};

/// Bytes [fileOffset, fileOffset + size) of the file belong at @p dest.
struct ByteRun {
  uint64_t fileOffset;
  uint64_t size;
  char *dest;
};

/// One read from the file: bytes [fileOffset, fileOffset + size), holding the runs in @p pieces.
struct SpanRead {
  uint64_t fileOffset;
  uint64_t size;
  std::vector<ByteRun> pieces;
};

/** The runs that hold elements [first, first + count) of @p column, written consecutively from @p dest.
 *
 * Stored chunks of an unfiltered dataset are the values themselves, so any range of elements within a chunk is one
 * contiguous piece of the file. Chunks that sit back to back in the file are merged into one run.
 */
MANTID_DATAHANDLING_DLL std::vector<ByteRun> planRuns(const ChunkedColumn &column, uint64_t first, uint64_t count,
                                                      char *dest);

/** Group runs into reads in file order.
 *
 * Runs are sorted by file offset, split to at most @p maxSpan bytes, and merged while the gap between them is at most
 * @p maxGap and the read stays within @p maxSpan. Each read is widened to whole blocks of @p block bytes (never past
 * @p fileSize), so a block shared by several runs is requested once.
 */
MANTID_DATAHANDLING_DLL std::vector<SpanRead> planSpans(std::vector<ByteRun> runs, uint64_t block, uint64_t maxSpan,
                                                        uint64_t maxGap, uint64_t fileSize);

/** Reads event columns of an event NeXus file straight from the file, without HDF5.
 *
 * HDF5 serialises every read behind one library-wide lock and keeps one request outstanding, and it finds chunks by
 * walking each dataset's index one small read at a time. In files written by the SNS data acquisition the index
 * nodes are scattered through the whole file, so on network storage that walk alone can take longer than reading
 * the data. This class parses the chunk indexes itself and reads the chunks directly on a pool of threads.
 *
 * Only the upper levels of each index are read up front, which gives every leaf node's address and the chunks it
 * covers. The leaves, nearly all of the index, are read when their chunks are first needed (locateChunks), so finding
 * the chunks overlaps reading the events instead of preceding it.
 *
 * Only what raw event files contain is handled: chunked, unfiltered, one-dimensional ``event_id`` (uint32) and
 * ``event_time_offset`` (float32) in native byte order, in files with a version 0 or 1 superblock and version 1
 * object headers. Any other column is left to HDF5.
 */
class MANTID_DATAHANDLING_DLL DirectEventReader {
public:
  DirectEventReader(const std::string &filename, H5::H5File &file, const std::vector<std::string> &bankEntryNames,
                    size_t numThreads = 32);

  /// The layout of the column at @p path (e.g. "/entry/bank1_events/event_id"), or nullptr if HDF5 must read it.
  /// Only the offsets of chunks passed to locateChunks are filled in.
  const ChunkedColumn *column(const std::string &path) const;
  /// Number of columns that will be read directly, out of the number looked at.
  size_t numDirectColumns() const { return m_columns.size(); }
  size_t numColumnsExamined() const { return m_numExamined; }

  /// Fill in the offsets of every chunk holding the given elements, reading any index leaves not yet read in one
  /// round of reads that go ahead of queued data reads. Safe to call from several threads.
  void locateChunks(const std::vector<ElementRange> &ranges) const;
  /// Roughly where the chunk holding @p element sits in the file, before it has been located: the address of the
  /// index leaf covering it, which the data acquisition writes next to the data. For ordering reads.
  uint64_t positionHint(const std::string &path, uint64_t element) const;

  /// Upper index levels read up front: nodes and seconds
  size_t numIndexNodes() const { return m_numIndexNodes; }
  double indexSeconds() const { return m_indexSeconds; }
  /// Index leaves read on demand so far, and the seconds spent waiting for them
  size_t numLeavesRead() const;
  double leafSeconds() const;

  /** Read the elements of each slab of the column at @p path consecutively into @p dest.
   *
   * @return false, having read nothing, if the column is not read directly or its element size is not
   * @p elementSize
   */
  bool read(const std::string &path, uint32_t elementSize, const std::vector<size_t> &offsets,
            const std::vector<size_t> &slabsizes, char *dest) const;

  ParallelFileReader &fileReader() const { return *m_reader; }
  /// The file system's block size, the unit it reads in.
  uint64_t blockSize() const { return m_blockSize; }

  /// An index leaf node: its address, and the chunks [firstChunk, endChunk) it holds
  struct Leaf {
    uint64_t address;
    uint64_t firstChunk;
    uint64_t endChunk;
  };
  /// A column's index leaves, in chunk order, and which have been read
  struct ColumnIndex {
    std::vector<Leaf> leaves;
    std::vector<char> read;
    uint32_t dimensionality{0};
    uint64_t nodeBytes{0};
  };

private:
  /// chunk offsets and which leaves are read are filled in by locateChunks, which is const so that shared readers can
  /// call it; m_indexMutex guards them
  mutable std::map<std::string, ChunkedColumn> m_columns;
  mutable std::map<std::string, ColumnIndex> m_indexes;
  mutable std::mutex m_indexMutex;
  mutable size_t m_numLeavesRead{0};
  mutable double m_leafSeconds{0.};
  uint32_t m_sizeOfOffsets{8};
  size_t m_numExamined{0};
  size_t m_numIndexNodes{0};
  double m_indexSeconds{0.};
  std::unique_ptr<ParallelFileReader> m_reader;
  uint64_t m_blockSize;
};

} // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
