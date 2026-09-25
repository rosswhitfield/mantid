// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include <cxxtest/TestSuite.h>

#include "MantidAPI/FileFinder.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/DirectEventReader.h"

#include <H5Cpp.h>
#include <Poco/TemporaryFile.h>
#include <algorithm>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>

using Mantid::DataHandling::AlignAndFocusPowderSlim::ByteRun;
using Mantid::DataHandling::AlignAndFocusPowderSlim::DirectEventReader;
using Mantid::DataHandling::AlignAndFocusPowderSlim::ElementRange;
using Mantid::DataHandling::AlignAndFocusPowderSlim::ParallelFileReader;
using Mantid::DataHandling::AlignAndFocusPowderSlim::planSpans;

namespace {
constexpr hsize_t CHUNK = 16;
// 2,000 chunks with a padded last one: enough for a B-tree more than one level deep
constexpr hsize_t NUM_EVENTS = CHUNK * 2000 - 5;
constexpr hsize_t PIECE = 100; // the columns are written alternately in pieces this long, so their chunks interleave

uint32_t detidValue(hsize_t i) { return static_cast<uint32_t>(i * 7 + 3); }
float tofValue(hsize_t i) { return static_cast<float>(i) * 0.5f + 0.25f; }

/// A growable one-dimensional chunked dataset, optionally compressed
H5::DataSet makeColumn(H5::Group &group, const std::string &name, const H5::PredType &type, bool deflate) {
  const hsize_t dims[1] = {0};
  const hsize_t maxdims[1] = {H5S_UNLIMITED};
  H5::DataSpace space(1, dims, maxdims);
  H5::DSetCreatPropList props;
  const hsize_t chunk[1] = {CHUNK};
  props.setChunk(1, chunk);
  if (deflate)
    props.setDeflate(4);
  return group.createDataSet(name, type, space, props);
}

template <typename T> void appendTo(H5::DataSet &dataset, const H5::PredType &type, const std::vector<T> &values) {
  hsize_t size[1];
  dataset.getSpace().getSimpleExtentDims(size);
  const hsize_t start[1] = {size[0]};
  const hsize_t count[1] = {values.size()};
  const hsize_t grown[1] = {size[0] + values.size()};
  dataset.extend(grown);
  auto filespace = dataset.getSpace();
  filespace.selectHyperslab(H5S_SELECT_SET, count, start);
  H5::DataSpace memspace(1, count);
  dataset.write(values.data(), type, memspace, filespace);
}

/// Write detector ids and times of flight to a bank alternately, the way the data acquisition interleaves them
void writeInterleaved(H5::Group &group, bool deflate) {
  auto detid = makeColumn(group, "event_id", H5::PredType::NATIVE_UINT32, deflate);
  auto tof = makeColumn(group, "event_time_offset", H5::PredType::NATIVE_FLOAT, deflate);
  for (hsize_t first = 0; first < NUM_EVENTS; first += PIECE) {
    const hsize_t last = std::min(first + PIECE, NUM_EVENTS);
    std::vector<uint32_t> ids;
    std::vector<float> tofs;
    for (hsize_t i = first; i < last; ++i) {
      ids.push_back(detidValue(i));
      tofs.push_back(tofValue(i));
    }
    appendTo(detid, H5::PredType::NATIVE_UINT32, ids);
    appendTo(tof, H5::PredType::NATIVE_FLOAT, tofs);
  }
}

/** A small event file:
 *  bank1: plain interleaved chunked columns (read directly)
 *  bank2: the same, compressed (left to HDF5)
 *  bank3: contiguous event_id (left to HDF5), chunked event_time_offset (read directly)
 *  bank4: float64 event_time_offset (left to HDF5), chunked event_id (read directly)
 *  bank5: empty chunked columns (read directly, nothing to locate)
 */
void writeTestFile(const std::string &filename) {
  H5::H5File file(filename, H5F_ACC_TRUNC);
  auto entry = file.createGroup("entry");
  auto bank1 = entry.createGroup("bank1_events");
  writeInterleaved(bank1, false);
  auto bank2 = entry.createGroup("bank2_events");
  writeInterleaved(bank2, true);

  auto bank3 = entry.createGroup("bank3_events");
  const hsize_t dims[1] = {100};
  std::vector<uint32_t> ids(100, 1);
  bank3.createDataSet("event_id", H5::PredType::NATIVE_UINT32, H5::DataSpace(1, dims))
      .write(ids.data(), H5::PredType::NATIVE_UINT32);
  auto tof3 = makeColumn(bank3, "event_time_offset", H5::PredType::NATIVE_FLOAT, false);
  appendTo(tof3, H5::PredType::NATIVE_FLOAT, std::vector<float>(100, 2.f));

  auto bank4 = entry.createGroup("bank4_events");
  auto id4 = makeColumn(bank4, "event_id", H5::PredType::NATIVE_UINT32, false);
  appendTo(id4, H5::PredType::NATIVE_UINT32, std::vector<uint32_t>(100, 3));
  auto tof4 = makeColumn(bank4, "event_time_offset", H5::PredType::NATIVE_DOUBLE, false);
  appendTo(tof4, H5::PredType::NATIVE_DOUBLE, std::vector<double>(100, 4.));

  auto bank5 = entry.createGroup("bank5_events");
  makeColumn(bank5, "event_id", H5::PredType::NATIVE_UINT32, false);
  makeColumn(bank5, "event_time_offset", H5::PredType::NATIVE_FLOAT, false);
}

herr_t collectChunk(const hsize_t *offset, unsigned, haddr_t address, hsize_t, void *data) {
  auto &found = *static_cast<std::vector<std::pair<hsize_t, uint64_t>> *>(data);
  found.emplace_back(offset[0], static_cast<uint64_t>(address));
  return H5_ITER_CONT;
}

/// Chunk offsets as HDF5 reports them, in element order
std::vector<uint64_t> chunkOffsetsFromHDF5(H5::H5File &file, const std::string &path) {
  auto dataset = file.openDataSet(path);
  std::vector<std::pair<hsize_t, uint64_t>> found;
  H5Dchunk_iter(dataset.getId(), H5P_DEFAULT, collectChunk, &found);
  std::sort(found.begin(), found.end());
  std::vector<uint64_t> offsets;
  for (const auto &chunk : found)
    offsets.push_back(chunk.second);
  return offsets;
}
} // namespace

class DirectEventReaderTest : public CxxTest::TestSuite {
public:
  static DirectEventReaderTest *createSuite() { return new DirectEventReaderTest(); }
  static void destroySuite(DirectEventReaderTest *suite) { delete suite; }

  DirectEventReaderTest() : m_banks{"bank1_events", "bank2_events", "bank3_events", "bank4_events", "bank5_events"} {
    writeTestFile(m_file.path());
  }

  void test_only_plain_chunked_columns_are_read_directly() {
    H5::H5File file(m_file.path(), H5F_ACC_RDONLY);
    DirectEventReader reader(m_file.path(), file, m_banks, 4);
    TS_ASSERT_EQUALS(reader.numColumnsExamined(), 10);
    TS_ASSERT(reader.column("/entry/bank1_events/event_id"));
    TS_ASSERT(reader.column("/entry/bank1_events/event_time_offset"));
    TS_ASSERT(!reader.column("/entry/bank2_events/event_id"));          // compressed
    TS_ASSERT(!reader.column("/entry/bank2_events/event_time_offset")); // compressed
    TS_ASSERT(!reader.column("/entry/bank3_events/event_id"));          // contiguous
    TS_ASSERT(reader.column("/entry/bank3_events/event_time_offset"));
    TS_ASSERT(reader.column("/entry/bank4_events/event_id"));
    TS_ASSERT(!reader.column("/entry/bank4_events/event_time_offset")); // float64
    TS_ASSERT(reader.column("/entry/bank5_events/event_id"));           // empty
    TS_ASSERT_EQUALS(reader.column("/entry/bank5_events/event_id")->numElements, 0);
    TS_ASSERT_EQUALS(reader.numDirectColumns(), 6);
  }

  void test_chunk_offsets_match_hdf5() {
    H5::H5File file(m_file.path(), H5F_ACC_RDONLY);
    DirectEventReader reader(m_file.path(), file, m_banks, 4);
    for (const std::string path : {"/entry/bank1_events/event_id", "/entry/bank1_events/event_time_offset"}) {
      reader.locateChunks({{path, 0, NUM_EVENTS}});
      const auto *column = reader.column(path);
      TS_ASSERT(column);
      TS_ASSERT_EQUALS(column->numElements, NUM_EVENTS);
      TS_ASSERT_EQUALS(column->chunkElements, CHUNK);
      TS_ASSERT_EQUALS(column->chunkOffsets, chunkOffsetsFromHDF5(file, path));
    }
  }

  void test_read_gives_the_values_written() {
    H5::H5File file(m_file.path(), H5F_ACC_RDONLY);
    DirectEventReader reader(m_file.path(), file, m_banks, 4);
    // whole columns, one partial chunk, a slab ending in the padded last chunk, and several slabs at once
    const std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> cases{
        {{0}, {NUM_EVENTS}}, {{5}, {7}}, {{NUM_EVENTS - 40}, {40}}, {{3, 100, 1000, 20000}, {50, 1, 999, 11000}}};
    for (const auto &[offsets, slabsizes] : cases) {
      const size_t total = std::accumulate(slabsizes.cbegin(), slabsizes.cend(), size_t{0});
      std::vector<uint32_t> ids(total);
      std::vector<float> tofs(total);
      TS_ASSERT(
          reader.read("/entry/bank1_events/event_id", 4, offsets, slabsizes, reinterpret_cast<char *>(ids.data())));
      TS_ASSERT(reader.read("/entry/bank1_events/event_time_offset", 4, offsets, slabsizes,
                            reinterpret_cast<char *>(tofs.data())));
      size_t position = 0;
      bool same = true;
      for (size_t slab = 0; slab < offsets.size(); ++slab) {
        for (size_t i = offsets[slab]; i < offsets[slab] + slabsizes[slab]; ++i, ++position)
          same = same && ids[position] == detidValue(i) && tofs[position] == tofValue(i);
      }
      TS_ASSERT(same);
    }
  }

  void test_index_leaves_are_read_in_the_background() {
    H5::H5File file(m_file.path(), H5F_ACC_RDONLY);
    DirectEventReader reader(m_file.path(), file, m_banks, 4);
    TS_ASSERT_LESS_THAN(0, reader.numIndexNodes()); // the upper levels, read up front
    const std::string idPath("/entry/bank1_events/event_id");
    const std::string tofPath("/entry/bank1_events/event_time_offset");
    const auto *ids = reader.column(idPath);
    const auto *tofs = reader.column(tofPath);

    // a hint is available before anything is located: a leaf's address or, once its leaf is read, the chunk's offset
    TS_ASSERT_LESS_THAN(0, reader.positionHint(idPath, 0));

    // locating waits for the leaves the range needs
    reader.locateChunks({{idPath, 5, 7}});
    TS_ASSERT_DIFFERS(ids->chunkOffsets[0], std::numeric_limits<uint64_t>::max());
    TS_ASSERT_EQUALS(reader.positionHint(idPath, 0), ids->chunkOffsets[0]);

    // once every chunk of bank1, the only columns with more than one leaf, is located, every queued leaf is read
    reader.locateChunks({{idPath, 0, NUM_EVENTS}, {tofPath, 0, NUM_EVENTS}});
    const auto unset = [](uint64_t offset) { return offset == std::numeric_limits<uint64_t>::max(); };
    TS_ASSERT(std::none_of(ids->chunkOffsets.cbegin(), ids->chunkOffsets.cend(), unset));
    TS_ASSERT(std::none_of(tofs->chunkOffsets.cbegin(), tofs->chunkOffsets.cend(), unset));
    TS_ASSERT_LESS_THAN(10, reader.numLeavesRead());
    TS_ASSERT_LESS_THAN(0., reader.leafBackgroundSeconds());
  }

  void test_destroying_the_reader_with_leaves_still_queued() {
    H5::H5File file(m_file.path(), H5F_ACC_RDONLY);
    for (int i = 0; i < 5; ++i) {
      DirectEventReader reader(m_file.path(), file, m_banks, 2); // destroyed before the leaves can all be read
    }
    TS_ASSERT(true);
  }

  void test_read_declines_what_it_does_not_handle() {
    H5::H5File file(m_file.path(), H5F_ACC_RDONLY);
    DirectEventReader reader(m_file.path(), file, m_banks, 4);
    std::vector<uint64_t> buffer(10);
    auto *dest = reinterpret_cast<char *>(buffer.data());
    TS_ASSERT(!reader.read("/entry/bank2_events/event_id", 4, {0}, {10}, dest)); // compressed
    TS_ASSERT(!reader.read("/entry/bank1_events/event_id", 8, {0}, {10}, dest)); // wrong element size
    TS_ASSERT(!reader.read("/entry/bank9_events/event_id", 4, {0}, {10}, dest)); // no such column
    TS_ASSERT_THROWS(reader.read("/entry/bank1_events/event_id", 4, {NUM_EVENTS}, {1}, dest),
                     const std::out_of_range &);
  }

  void test_chunk_offsets_match_hdf5_for_a_real_file() {
    const auto filename = Mantid::API::FileFinder::Instance().getFullPath("VULCAN_218062.nxs.h5");
    TS_ASSERT(!filename.empty());
    H5::H5File file(filename, H5F_ACC_RDONLY);
    const std::vector<std::string> banks{"bank1_events", "bank2_events", "bank3_events",
                                         "bank4_events", "bank5_events", "bank6_events"};
    DirectEventReader reader(filename, file, banks);
    TS_ASSERT_EQUALS(reader.numDirectColumns(), 12);
    for (const auto &bank : banks) {
      for (const std::string name : {"event_id", "event_time_offset"}) {
        const auto path = "/entry/" + bank + "/" + name;
        const auto *column = reader.column(path);
        TS_ASSERT(column);
        if (column) {
          reader.locateChunks({{path, 0, column->numElements}});
          TS_ASSERT_EQUALS(column->chunkOffsets, chunkOffsetsFromHDF5(file, path));
        }
      }
    }
  }

  void test_spans_merge_through_small_gaps_and_align_to_blocks() {
    std::vector<char> dest(100);
    char *base = dest.data();
    // two runs 10 bytes apart, then one far away
    const std::vector<ByteRun> runs{{1030, 20, base}, {1060, 5, base + 20}, {9000, 10, base + 30}};
    const auto spans = planSpans(runs, 512, 4096, 512, 100000);
    TS_ASSERT_EQUALS(spans.size(), 2);
    TS_ASSERT_EQUALS(spans[0].fileOffset, 1024); // widened to whole 512-byte blocks
    TS_ASSERT_EQUALS(spans[0].size, 512);
    TS_ASSERT_EQUALS(spans[0].pieces.size(), 2);
    TS_ASSERT_EQUALS(spans[1].fileOffset, 8704);
    TS_ASSERT_EQUALS(spans[1].size, 512);
    // no gap tolerance: the first two runs are read separately
    TS_ASSERT_EQUALS(planSpans(runs, 1, 4096, 0, 100000).size(), 3);
  }

  void test_spans_split_long_runs_and_stop_at_the_end_of_the_file() {
    std::vector<char> dest(10000);
    const std::vector<ByteRun> runs{{100, 10000, dest.data()}};
    // the file ends at 10150, inside the block the last run ends in
    const auto spans = planSpans(runs, 1000, 4000, 1000, 10150);
    uint64_t covered = 0;
    for (const auto &span : spans) {
      TS_ASSERT_LESS_THAN_EQUALS(span.size, 4000 + 2 * 1000);
      TS_ASSERT_LESS_THAN_EQUALS(span.fileOffset + span.size, 10150);
      for (const auto &piece : span.pieces) {
        TS_ASSERT_LESS_THAN_EQUALS(span.fileOffset, piece.fileOffset);
        TS_ASSERT_LESS_THAN_EQUALS(piece.fileOffset + piece.size, span.fileOffset + span.size);
        // each piece keeps its place in the destination
        TS_ASSERT_EQUALS(piece.dest - dest.data(), static_cast<std::ptrdiff_t>(piece.fileOffset - 100));
        covered += piece.size;
      }
    }
    TS_ASSERT_EQUALS(covered, 10000);
  }

  void test_pool_reports_the_first_failure() {
    ParallelFileReader pool(m_file.path(), 3);
    std::vector<ParallelFileReader::Task> tasks;
    for (int i = 0; i < 10; ++i)
      tasks.emplace_back([i](std::istream &) {
        if (i == 4)
          throw std::runtime_error("task 4 failed");
      });
    TS_ASSERT_THROWS(pool.run(std::move(tasks)), const std::runtime_error &);
    // and keeps working afterwards
    std::vector<char> head(8);
    pool.run({[&head](std::istream &stream) {
      Mantid::DataHandling::AlignAndFocusPowderSlim::readAt(stream, 0, 8, head.data());
    }});
    TS_ASSERT_EQUALS(std::memcmp(head.data(), "\x89HDF\r\n\x1a\n", 8), 0);
  }

private:
  Poco::TemporaryFile m_file;
  const std::vector<std::string> m_banks;
};
