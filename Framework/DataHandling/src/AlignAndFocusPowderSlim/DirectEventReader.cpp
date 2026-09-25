// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#include "MantidDataHandling/AlignAndFocusPowderSlim/DirectEventReader.h"
#include "MantidKernel/Logger.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>

#ifndef _WIN32
#include <sys/statvfs.h>
#endif

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

namespace {
auto g_log = Kernel::Logger("DirectEventReader");

// ---------------------------------------------------------------------------
// HDF5 file format constants (HDF5 File Format Specification, version 3.0)

const char HDF5_SIGNATURE[8] = {'\x89', 'H', 'D', 'F', '\r', '\n', '\x1a', '\n'};
constexpr uint16_t LAYOUT_MESSAGE = 0x0008;
constexpr uint16_t CONTINUATION_MESSAGE = 0x0010;
constexpr uint8_t LAYOUT_CHUNKED = 2;
constexpr uint8_t BTREE_RAW_DATA_NODE = 1;
/// Indexed-storage B-tree K when the superblock does not record one (version 0)
constexpr uint32_t DEFAULT_ISTORE_K = 32;
constexpr uint64_t DEFAULT_BLOCK_SIZE = 4 * 1024 * 1024;
/// Largest piece a single read of chunk data is split into
constexpr uint64_t MAX_DIRECT_READ = 16 * 1024 * 1024;

/// The file format uses something this reader does not handle; the column is left to HDF5.
class Unsupported : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

uint64_t readLittleEndian(const char *data, size_t numBytes) {
  uint64_t value = 0;
  for (size_t i = 0; i < numBytes; ++i)
    value |= static_cast<uint64_t>(static_cast<unsigned char>(data[i])) << (8 * i);
  return value;
}

/// Read up to @p size bytes at @p offset; fewer at the end of the file.
std::vector<char> readUpTo(std::istream &stream, uint64_t offset, uint64_t size) {
  std::vector<char> buffer(size);
  stream.clear();
  stream.seekg(static_cast<std::streamoff>(offset));
  stream.read(buffer.data(), static_cast<std::streamsize>(size));
  buffer.resize(static_cast<size_t>(stream.gcount()));
  stream.clear();
  return buffer;
}

struct Superblock {
  uint32_t sizeOfOffsets;
  uint32_t istoreK;
};

Superblock readSuperblock(std::istream &stream) {
  const auto head = readUpTo(stream, 0, 64);
  if (head.size() < 40 || std::memcmp(head.data(), HDF5_SIGNATURE, 8) != 0)
    throw Unsupported("no superblock at the start of the file (user block?)");
  const auto version = static_cast<uint8_t>(head[8]);
  if (version > 1)
    throw Unsupported("superblock version " + std::to_string(version));
  Superblock superblock{static_cast<uint8_t>(head[13]), DEFAULT_ISTORE_K};
  if (superblock.sizeOfOffsets != 4 && superblock.sizeOfOffsets != 8)
    throw Unsupported("size of offsets " + std::to_string(superblock.sizeOfOffsets));
  // version 1 records the indexed storage K after the group Ks and file consistency flags
  size_t baseAddressPosition = 24;
  if (version == 1) {
    superblock.istoreK = static_cast<uint32_t>(readLittleEndian(head.data() + 24, 2));
    baseAddressPosition = 28;
  }
  if (readLittleEndian(head.data() + baseAddressPosition, superblock.sizeOfOffsets) != 0)
    throw Unsupported("non-zero base address");
  return superblock;
}

/// The root address of a chunked dataset's version-1 B-tree and its dimensionality (rank + 1), from the layout
/// message of its version-1 object header.
std::pair<uint64_t, uint32_t> btreeRoot(std::istream &stream, uint64_t headerAddress, const Superblock &superblock) {
  const auto head = readUpTo(stream, headerAddress, 16);
  if (head.size() < 16 || head[0] != 1)
    throw Unsupported("object header is not version 1");
  const auto numMessages = readLittleEndian(head.data() + 2, 2);
  std::deque<std::pair<uint64_t, uint64_t>> blocks{{headerAddress + 16, readLittleEndian(head.data() + 8, 4)}};
  uint64_t seen = 0;
  const auto offsetSize = superblock.sizeOfOffsets;
  while (!blocks.empty() && seen < numMessages) {
    const auto [start, length] = blocks.front();
    blocks.pop_front();
    const auto data = readUpTo(stream, start, length);
    size_t position = 0;
    while (position + 8 <= data.size() && seen < numMessages) {
      const auto type = static_cast<uint16_t>(readLittleEndian(data.data() + position, 2));
      const auto size = readLittleEndian(data.data() + position + 2, 2);
      const char *body = data.data() + position + 8;
      if (position + 8 + size > data.size())
        throw Unsupported("truncated object header message");
      ++seen;
      if (type == CONTINUATION_MESSAGE) {
        blocks.emplace_back(readLittleEndian(body, offsetSize), readLittleEndian(body + offsetSize, offsetSize));
      } else if (type == LAYOUT_MESSAGE) {
        const auto version = static_cast<uint8_t>(body[0]);
        if (version == 1 || version == 2) {
          if (static_cast<uint8_t>(body[2]) != LAYOUT_CHUNKED)
            throw Unsupported("not chunked");
          return {readLittleEndian(body + 8, offsetSize), static_cast<uint8_t>(body[1])};
        }
        if (version == 3) {
          if (static_cast<uint8_t>(body[1]) != LAYOUT_CHUNKED)
            throw Unsupported("not chunked");
          return {readLittleEndian(body + 3, offsetSize), static_cast<uint8_t>(body[2])};
        }
        throw Unsupported("layout message version " + std::to_string(version) + " (newer chunk index types)");
      }
      position += 8 + size;
    }
  }
  throw Unsupported("no layout message");
}

/// A B-tree node of a raw-data chunk index, parsed
struct Node {
  uint8_t level{0};
  /// element offset of each key; one more than there are children, the last being an upper bound
  std::vector<uint64_t> keys;
  /// stored size of each child chunk (meaningful in leaves)
  std::vector<uint64_t> storedSizes;
  std::vector<uint64_t> children;
};

std::optional<Node> parseNode(const std::vector<char> &data, uint32_t dimensionality, uint64_t offsetSize) {
  if (data.size() < 8 + 2 * offsetSize || std::memcmp(data.data(), "TREE", 4) != 0 ||
      static_cast<uint8_t>(data[4]) != BTREE_RAW_DATA_NODE)
    return std::nullopt;
  Node node;
  node.level = static_cast<uint8_t>(data[5]);
  const auto entries = readLittleEndian(data.data() + 6, 2);
  const uint64_t keySize = 8 + 8 * static_cast<uint64_t>(dimensionality);
  if (8 + 2 * offsetSize + (entries + 1) * keySize + entries * offsetSize > data.size())
    return std::nullopt;
  size_t position = 8 + 2 * offsetSize; // skip the sibling addresses
  for (uint64_t entry = 0; entry <= entries; ++entry) {
    // key: stored chunk size (4 bytes), filter mask (4 bytes), then the chunk's offset in each dimension
    node.storedSizes.push_back(readLittleEndian(data.data() + position, 4));
    node.keys.push_back(readLittleEndian(data.data() + position + 8, 8));
    position += keySize;
    if (entry < entries) {
      node.children.push_back(readLittleEndian(data.data() + position, offsetSize));
      position += offsetSize;
    }
  }
  return node;
}

/// Record a leaf's chunk offsets; false if it is not the plain unfiltered leaf expected or leaves gaps in
/// [firstChunk, endChunk).
bool applyLeaf(const Node &leaf, ChunkedColumn &column, uint64_t firstChunk, uint64_t endChunk) {
  if (leaf.level != 0)
    return false;
  uint64_t found = 0;
  for (size_t entry = 0; entry < leaf.children.size(); ++entry) {
    const auto elementOffset = leaf.keys[entry];
    const auto chunk = elementOffset / column.chunkElements;
    if (elementOffset % column.chunkElements != 0 || chunk < firstChunk || chunk >= endChunk ||
        leaf.storedSizes[entry] != column.chunkElements * column.elementSize)
      return false;
    if (column.chunkOffsets[chunk] == std::numeric_limits<uint64_t>::max())
      ++found;
    column.chunkOffsets[chunk] = leaf.children[entry];
  }
  return found == endChunk - firstChunk; // otherwise some chunks were never written
}

/// A column whose chunk index is being read.
struct IndexJob {
  std::string path;
  ChunkedColumn column;
  uint32_t dimensionality{0};
  uint64_t nodeBytes{0};
  /// (address, first chunk, end chunk) of each leaf; a leaf read already has its chunks in column
  std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> leaves;
  std::vector<char> leafRead;
  bool failed{false};
};

/** Read the levels of every job's B-tree above the leaves, level by level: every node of a level, across all the
 * trees, is fetched at once on the pool. The leaves' addresses and chunk ranges are recorded; a tree whose root is
 * itself a leaf is read completely. Jobs whose trees cannot be read are marked failed.
 * @return the number of nodes read
 */
size_t readUpperLevels(ParallelFileReader &reader, std::vector<IndexJob> &jobs, std::vector<uint64_t> roots,
                       const Superblock &superblock) {
  size_t numNodes = 0;
  const uint64_t offsetSize = superblock.sizeOfOffsets;
  for (auto &job : jobs) {
    // a node holds at most 2K entries; keys are the chunk size, filter mask and one 8-byte offset per dimension
    const uint64_t keySize = 8 + 8 * static_cast<uint64_t>(job.dimensionality);
    job.nodeBytes = 8 + 2 * offsetSize + (2 * superblock.istoreK + 1) * keySize + 2 * superblock.istoreK * offsetSize;
  }

  std::vector<std::pair<size_t, uint64_t>> frontier; // (job, node address)
  for (size_t i = 0; i < jobs.size(); ++i)
    if (!jobs[i].failed)
      frontier.emplace_back(i, roots[i]);

  while (!frontier.empty()) {
    std::vector<std::vector<char>> nodes(frontier.size());
    std::vector<ParallelFileReader::Task> tasks;
    tasks.reserve(frontier.size());
    for (size_t n = 0; n < frontier.size(); ++n) {
      tasks.emplace_back([&, n](std::istream &stream) {
        nodes[n] = readUpTo(stream, frontier[n].second, jobs[frontier[n].first].nodeBytes);
      });
    }
    reader.run(std::move(tasks));
    numNodes += frontier.size();

    std::vector<std::pair<size_t, uint64_t>> next;
    for (size_t n = 0; n < frontier.size(); ++n) {
      auto &job = jobs[frontier[n].first];
      if (job.failed)
        continue;
      const auto node = parseNode(nodes[n], job.dimensionality, offsetSize);
      if (!node) {
        job.failed = true;
        continue;
      }
      const uint64_t numChunks = job.column.chunkOffsets.size();
      if (node->level == 0) { // the root is a leaf: read it now
        job.leaves.emplace_back(frontier[n].second, 0, numChunks);
        job.leafRead.push_back(1);
        job.failed = !applyLeaf(*node, job.column, 0, numChunks);
      } else if (node->level == 1) { // the children are leaves: record where they are and what they hold
        // each key bounds its child from below; the ends are set below from where the next leaf starts, since the
        // tree's last key is the last chunk itself rather than a bound past it
        for (size_t i = 0; i < node->children.size(); ++i) {
          job.leaves.emplace_back(node->children[i], node->keys[i] / job.column.chunkElements, numChunks);
          job.leafRead.push_back(0);
        }
      } else {
        for (const auto child : node->children)
          next.emplace_back(frontier[n].first, child);
      }
    }
    frontier = std::move(next);
  }

  // order the leaves by the chunks they hold; each runs up to where the next starts, the last to the end. Reading a
  // leaf checks that its chunks lie in that range and that none is missing.
  for (auto &job : jobs) {
    if (job.failed || job.leaves.empty())
      continue;
    std::vector<size_t> order(job.leaves.size());
    std::iota(order.begin(), order.end(), size_t{0});
    std::sort(order.begin(), order.end(),
              [&job](size_t a, size_t b) { return std::get<1>(job.leaves[a]) < std::get<1>(job.leaves[b]); });
    decltype(job.leaves) leaves;
    std::vector<char> leafRead;
    for (const auto i : order) {
      leaves.push_back(job.leaves[i]);
      leafRead.push_back(job.leafRead[i]);
    }
    const uint64_t numChunks = job.column.chunkOffsets.size();
    if (std::get<1>(leaves.front()) != 0)
      job.failed = true;
    for (size_t i = 0; i < leaves.size() && !job.failed; ++i) {
      const uint64_t end = i + 1 < leaves.size() ? std::get<1>(leaves[i + 1]) : numChunks;
      if (end <= std::get<1>(leaves[i]))
        job.failed = true; // two leaves claim the same first chunk
      std::get<2>(leaves[i]) = end;
    }
    job.leaves = std::move(leaves);
    job.leafRead = std::move(leafRead);
  }
  return numNodes;
}

uint64_t fileSystemBlockSize(const std::string &filename) {
#ifndef _WIN32
  struct statvfs info{};
  if (statvfs(filename.c_str(), &info) == 0 && info.f_bsize >= 4096 && info.f_bsize <= 64 * 1024 * 1024)
    return static_cast<uint64_t>(info.f_bsize);
#else
  (void)filename;
#endif
  return DEFAULT_BLOCK_SIZE;
}

/// The chunk layout of a dataset if it can be read directly: chunked, unfiltered, one-dimensional, of
/// @p expectedType. The chunk offsets are left to be filled in.
std::optional<IndexJob> describeColumn(H5::DataSet &dataset, const std::string &path, const H5::PredType &expectedType,
                                       uint64_t &headerAddress) {
  const hid_t id = dataset.getId();
  const hid_t dcpl = H5Dget_create_plist(id);
  hsize_t chunkDims[1] = {0};
  const bool chunked = H5Pget_layout(dcpl) == H5D_CHUNKED && H5Pget_nfilters(dcpl) == 0 &&
                       H5Pget_chunk(dcpl, 1, chunkDims) == 1 && chunkDims[0] > 0;
  H5Pclose(dcpl);
  if (!chunked)
    return std::nullopt;

  const hid_t space = H5Dget_space(id);
  hsize_t dims[1] = {0};
  const bool oneDimensional =
      H5Sget_simple_extent_ndims(space) == 1 && H5Sget_simple_extent_dims(space, dims, nullptr) == 1;
  H5Sclose(space);
  const hid_t type = H5Dget_type(id);
  const bool rightType = H5Tequal(type, expectedType.getId()) > 0;
  H5Tclose(type);
  if (!oneDimensional || !rightType)
    return std::nullopt;

  // basic fields only: the default H5Oget_info fields include the metadata size, which walks the whole chunk index
  H5O_info2_t info;
  haddr_t address = HADDR_UNDEF;
  if (H5Oget_info3(id, &info, H5O_INFO_BASIC) < 0 || H5VLnative_token_to_addr(id, info.token, &address) < 0)
    return std::nullopt;
  headerAddress = static_cast<uint64_t>(address);

  IndexJob job;
  job.path = path;
  job.column.numElements = dims[0];
  job.column.chunkElements = chunkDims[0];
  job.column.elementSize = static_cast<uint32_t>(expectedType.getSize());
  job.column.chunkOffsets.assign((dims[0] + chunkDims[0] - 1) / chunkDims[0], std::numeric_limits<uint64_t>::max());
  return job;
}
} // namespace

// ---------------------------------------------------------------------------
// ParallelFileReader

ParallelFileReader::ParallelFileReader(const std::string &filename, size_t numThreads)
    : m_filename(filename), m_fileSize(std::filesystem::file_size(filename)) {
  numThreads = std::max<size_t>(numThreads, 1);
  m_threads.reserve(numThreads);
  for (size_t i = 0; i < numThreads; ++i)
    m_threads.emplace_back([this] { work(); });
}

ParallelFileReader::~ParallelFileReader() {
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_stop = true;
  }
  m_wake.notify_all();
  for (auto &thread : m_threads)
    thread.join();
}

void ParallelFileReader::work() {
  std::ifstream stream(m_filename, std::ios::binary);
  while (true) {
    std::function<void(std::istream &)> task;
    {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_wake.wait(lock, [this] { return m_stop || !m_queue.empty(); });
      if (m_queue.empty())
        return; // stopping, and nothing left to do
      task = std::move(m_queue.front());
      m_queue.pop_front();
    }
    task(stream);
  }
}

std::future<void> ParallelFileReader::start(std::vector<Task> tasks, bool urgent) {
  struct Batch {
    std::promise<void> done;
    std::atomic<size_t> remaining;
    std::mutex errorMutex;
    std::exception_ptr error;
  };
  auto batch = std::make_shared<Batch>();
  auto future = batch->done.get_future();
  if (tasks.empty()) {
    batch->done.set_value();
    return future;
  }
  batch->remaining = tasks.size();
  {
    std::lock_guard<std::mutex> lock(m_mutex);
    // urgent tasks go to the front, in reverse so that they still run in the order given
    if (urgent)
      std::reverse(tasks.begin(), tasks.end());
    for (auto &task : tasks) {
      auto wrapped = [batch, task = std::move(task)](std::istream &stream) {
        try {
          task(stream);
        } catch (...) {
          std::lock_guard<std::mutex> errorLock(batch->errorMutex);
          if (!batch->error)
            batch->error = std::current_exception();
        }
        if (--batch->remaining == 0) {
          if (batch->error)
            batch->done.set_exception(batch->error);
          else
            batch->done.set_value();
        }
      };
      if (urgent)
        m_queue.emplace_front(std::move(wrapped));
      else
        m_queue.emplace_back(std::move(wrapped));
    }
  }
  m_wake.notify_all();
  return future;
}

void ParallelFileReader::run(std::vector<Task> tasks, bool urgent) { start(std::move(tasks), urgent).get(); }

void readAt(std::istream &stream, uint64_t offset, uint64_t size, char *dest) {
  stream.clear();
  stream.seekg(static_cast<std::streamoff>(offset));
  stream.read(dest, static_cast<std::streamsize>(size));
  if (static_cast<uint64_t>(stream.gcount()) != size) {
    stream.clear();
    throw std::runtime_error("Short read of " + std::to_string(size) + " bytes at offset " + std::to_string(offset));
  }
}

// ---------------------------------------------------------------------------
// planning

std::vector<ByteRun> planRuns(const ChunkedColumn &column, uint64_t first, uint64_t count, char *dest) {
  std::vector<ByteRun> runs;
  if (count == 0)
    return runs;
  if (first + count > column.numElements)
    throw std::out_of_range("Elements [" + std::to_string(first) + ", " + std::to_string(first + count) +
                            ") are beyond the end of a column of " + std::to_string(column.numElements));
  const uint64_t size = column.elementSize;
  const uint64_t chunkElements = column.chunkElements;
  const uint64_t last = first + count;
  for (uint64_t chunk = first / chunkElements; chunk * chunkElements < last; ++chunk) {
    const uint64_t chunkFirst = chunk * chunkElements;
    const uint64_t lo = std::max(first, chunkFirst);
    const uint64_t hi = std::min(last, chunkFirst + chunkElements);
    const uint64_t fileOffset = column.chunkOffsets[chunk] + (lo - chunkFirst) * size;
    const uint64_t bytes = (hi - lo) * size;
    char *target = dest + (lo - first) * size;
    if (!runs.empty() && runs.back().fileOffset + runs.back().size == fileOffset) {
      runs.back().size += bytes; // back to back in the file, and consecutive in the destination
    } else {
      runs.push_back({fileOffset, bytes, target});
    }
  }
  return runs;
}

std::vector<SpanRead> planSpans(std::vector<ByteRun> runs, uint64_t block, uint64_t maxSpan, uint64_t maxGap,
                                uint64_t fileSize) {
  block = std::max<uint64_t>(block, 1);
  maxSpan = std::max(maxSpan, block);
  // split long runs so that no read, and no read's staging buffer, exceeds maxSpan
  std::vector<ByteRun> pieces;
  pieces.reserve(runs.size());
  for (const auto &run : runs) {
    for (uint64_t done = 0; done < run.size; done += maxSpan)
      pieces.push_back({run.fileOffset + done, std::min(maxSpan, run.size - done), run.dest + done});
  }
  std::sort(pieces.begin(), pieces.end(),
            [](const ByteRun &left, const ByteRun &right) { return left.fileOffset < right.fileOffset; });

  std::vector<SpanRead> spans;
  uint64_t spanEnd = 0;
  for (const auto &piece : pieces) {
    const uint64_t pieceEnd = piece.fileOffset + piece.size;
    if (!spans.empty() && piece.fileOffset <= spanEnd + maxGap && pieceEnd - spans.back().fileOffset <= maxSpan) {
      spans.back().pieces.push_back(piece);
      spanEnd = std::max(spanEnd, pieceEnd);
      spans.back().size = spanEnd - spans.back().fileOffset;
    } else {
      spans.push_back({piece.fileOffset, piece.size, {piece}});
      spanEnd = pieceEnd;
    }
  }
  for (auto &span : spans) {
    const uint64_t start = span.fileOffset / block * block;
    const uint64_t end = std::min((span.fileOffset + span.size + block - 1) / block * block, fileSize);
    span.fileOffset = start;
    span.size = end - start;
  }
  return spans;
}

// ---------------------------------------------------------------------------
// DirectEventReader

DirectEventReader::DirectEventReader(const std::string &filename, H5::H5File &file,
                                     const std::vector<std::string> &bankEntryNames, size_t numThreads)
    : m_reader(std::make_unique<ParallelFileReader>(filename, numThreads)), m_blockSize(fileSystemBlockSize(filename)) {
  std::vector<IndexJob> jobs;
  std::vector<uint64_t> roots;
  std::vector<uint64_t> headers;
  const std::vector<std::pair<std::string, const H5::PredType *>> columns{
      {"event_id", &H5::PredType::NATIVE_UINT32}, {"event_time_offset", &H5::PredType::NATIVE_FLOAT}};
  for (const auto &bank : bankEntryNames) {
    if (bank.empty())
      continue;
    for (const auto &[name, type] : columns) {
      const std::string path = "/entry/" + bank + "/" + name;
      ++m_numExamined;
      try {
        auto dataset = file.openDataSet(path);
        uint64_t header = 0;
        auto job = describeColumn(dataset, path, *type, header);
        if (!job) {
          g_log.information() << path << " is not an unfiltered chunked column of the expected type; HDF5 reads it\n";
          continue;
        }
        if (job->column.numElements == 0) {
          m_columns.emplace(path, std::move(job->column)); // empty: nothing to locate
          continue;
        }
        jobs.push_back(std::move(*job));
        headers.push_back(header);
      } catch (const H5::Exception &) {
        g_log.information() << "could not open " << path << "; HDF5 reads it\n";
      }
    }
  }
  if (jobs.empty())
    return;

  const auto start = std::chrono::steady_clock::now();
  try {
    std::ifstream stream(filename, std::ios::binary);
    const auto superblock = readSuperblock(stream);
    for (size_t i = 0; i < jobs.size(); ++i) {
      try {
        const auto [root, dimensionality] = btreeRoot(stream, headers[i], superblock);
        jobs[i].dimensionality = dimensionality;
        roots.push_back(root);
      } catch (const Unsupported &error) {
        g_log.information() << jobs[i].path << ": " << error.what() << "; HDF5 reads it\n";
        jobs[i].failed = true;
        roots.push_back(0);
      }
    }
    m_sizeOfOffsets = superblock.sizeOfOffsets;
    m_numIndexNodes = readUpperLevels(*m_reader, jobs, std::move(roots), superblock);
    m_indexSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  } catch (const Unsupported &error) {
    g_log.information() << filename << ": " << error.what() << "; HDF5 reads every column\n";
    return;
  }
  for (auto &job : jobs) {
    if (job.failed) {
      g_log.information() << job.path << ": chunk index not understood; HDF5 reads it\n";
      continue;
    }
    ColumnIndex index;
    index.dimensionality = job.dimensionality;
    index.nodeBytes = job.nodeBytes;
    for (const auto &[address, first, end] : job.leaves)
      index.leaves.push_back({address, first, end});
    index.read = std::move(job.leafRead);
    m_indexes.emplace(job.path, std::move(index));
    m_columns.emplace(job.path, std::move(job.column));
  }
}

void DirectEventReader::locateChunks(const std::vector<ElementRange> &ranges) const {
  struct Wanted {
    std::string path;
    size_t leaf;
  };
  std::lock_guard<std::mutex> lock(m_indexMutex);
  std::vector<Wanted> wanted;
  std::set<std::pair<std::string, size_t>> seen;
  for (const auto &range : ranges) {
    if (range.count == 0)
      continue;
    const auto index = m_indexes.find(range.path);
    const auto column = m_columns.find(range.path);
    if (index == m_indexes.end() || column == m_columns.end())
      continue; // not read directly, or empty
    const auto &leaves = index->second.leaves;
    const uint64_t firstChunk = range.first / column->second.chunkElements;
    const uint64_t lastChunk = (range.first + range.count - 1) / column->second.chunkElements;
    // the first leaf holding firstChunk, then every leaf up to lastChunk
    auto leaf = std::upper_bound(leaves.begin(), leaves.end(), firstChunk,
                                 [](uint64_t chunk, const Leaf &candidate) { return chunk < candidate.firstChunk; });
    if (leaf != leaves.begin())
      --leaf;
    for (; leaf != leaves.end() && leaf->firstChunk <= lastChunk; ++leaf) {
      const auto position = static_cast<size_t>(leaf - leaves.begin());
      if (!index->second.read[position] && seen.emplace(range.path, position).second)
        wanted.push_back({range.path, position});
    }
  }
  if (wanted.empty())
    return;

  const auto start = std::chrono::steady_clock::now();
  std::vector<std::vector<char>> nodes(wanted.size());
  std::vector<ParallelFileReader::Task> tasks;
  tasks.reserve(wanted.size());
  for (size_t n = 0; n < wanted.size(); ++n) {
    const auto &index = m_indexes.at(wanted[n].path);
    const auto address = index.leaves[wanted[n].leaf].address;
    const auto bytes = index.nodeBytes;
    tasks.emplace_back(
        [&nodes, n, address, bytes](std::istream &stream) { nodes[n] = readUpTo(stream, address, bytes); });
  }
  // ahead of any queued data reads, which may be waiting on these very chunks
  m_reader->run(std::move(tasks), true);

  for (size_t n = 0; n < wanted.size(); ++n) {
    auto &index = m_indexes.at(wanted[n].path);
    const auto &leaf = index.leaves[wanted[n].leaf];
    const auto node = parseNode(nodes[n], index.dimensionality, m_sizeOfOffsets);
    if (!node || !applyLeaf(*node, m_columns.at(wanted[n].path), leaf.firstChunk, leaf.endChunk))
      throw std::runtime_error("Cannot read the chunk index of " + wanted[n].path + " at offset " +
                               std::to_string(leaf.address));
    index.read[wanted[n].leaf] = 1;
  }
  m_numLeavesRead += wanted.size();
  m_leafSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

uint64_t DirectEventReader::positionHint(const std::string &path, uint64_t element) const {
  std::lock_guard<std::mutex> lock(m_indexMutex);
  const auto index = m_indexes.find(path);
  const auto column = m_columns.find(path);
  if (index == m_indexes.end() || column == m_columns.end() || column->second.chunkOffsets.empty())
    return 0;
  const uint64_t chunk =
      std::min<uint64_t>(element / column->second.chunkElements, column->second.chunkOffsets.size() - 1);
  const auto &leaves = index->second.leaves;
  auto leaf = std::upper_bound(leaves.begin(), leaves.end(), chunk,
                               [](uint64_t value, const Leaf &candidate) { return value < candidate.firstChunk; });
  if (leaf != leaves.begin())
    --leaf;
  const auto position = static_cast<size_t>(leaf - leaves.begin());
  return index->second.read[position] ? column->second.chunkOffsets[chunk] : leaf->address;
}

size_t DirectEventReader::numLeavesRead() const {
  std::lock_guard<std::mutex> lock(m_indexMutex);
  return m_numLeavesRead;
}

double DirectEventReader::leafSeconds() const {
  std::lock_guard<std::mutex> lock(m_indexMutex);
  return m_leafSeconds;
}

const ChunkedColumn *DirectEventReader::column(const std::string &path) const {
  const auto found = m_columns.find(path);
  return found == m_columns.end() ? nullptr : &found->second;
}

bool DirectEventReader::read(const std::string &path, uint32_t elementSize, const std::vector<size_t> &offsets,
                             const std::vector<size_t> &slabsizes, char *dest) const {
  const auto *layout = column(path);
  if (layout == nullptr || layout->elementSize != elementSize)
    return false;
  std::vector<ElementRange> ranges;
  for (size_t i = 0; i < offsets.size(); ++i)
    ranges.push_back({path, offsets[i], slabsizes[i]});
  locateChunks(ranges);
  std::vector<ParallelFileReader::Task> tasks;
  char *target = dest;
  for (size_t i = 0; i < offsets.size(); ++i) {
    for (const auto &run : planRuns(*layout, offsets[i], slabsizes[i], target)) {
      // split long runs so several threads share them
      for (uint64_t done = 0; done < run.size; done += MAX_DIRECT_READ) {
        const uint64_t size = std::min(MAX_DIRECT_READ, run.size - done);
        tasks.emplace_back([fileOffset = run.fileOffset + done, size, out = run.dest + done](std::istream &stream) {
          readAt(stream, fileOffset, size, out);
        });
      }
    }
    target += static_cast<uint64_t>(slabsizes[i]) * elementSize;
  }
  m_reader->run(std::move(tasks));
  return true;
}

} // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
