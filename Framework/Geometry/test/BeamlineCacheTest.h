// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include <cxxtest/TestSuite.h>

#include "MantidFrameworkTestHelpers/ComponentCreationHelper.h"
#include "MantidGeometry/Instrument/BeamlineCache.h"
#include "MantidGeometry/Instrument/ComponentInfo.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"
#include "MantidGeometry/Instrument/InstrumentVisitor.h"
#include "MantidKernel/V3D.h"

#include <filesystem>
#include <fstream>

using namespace Mantid::Geometry;
using Mantid::detid_t;
using Mantid::Kernel::V3D;

namespace {

/// A cache file path in the system temporary directory, removed by the
/// destructor so that a failing assertion cannot leave a file behind.
class TemporaryCacheFile {
public:
  explicit TemporaryCacheFile(const std::string &name)
      : m_path((std::filesystem::temp_directory_path() / name).string()) {
    remove();
  }
  ~TemporaryCacheFile() { remove(); }
  TemporaryCacheFile(const TemporaryCacheFile &) = delete;
  TemporaryCacheFile &operator=(const TemporaryCacheFile &) = delete;

  const std::string &path() const { return m_path; }
  bool exists() const { return std::filesystem::exists(m_path); }

private:
  void remove() {
    std::error_code ec;
    std::filesystem::remove(m_path, ec);
    for (const auto &entry : std::filesystem::directory_iterator(std::filesystem::temp_directory_path(), ec))
      if (entry.path().string().rfind(m_path + ".tmp", 0) == 0)
        std::filesystem::remove(entry.path(), ec);
  }
  std::string m_path;
};

BeamlineCacheData makeData(const size_t numDetectors, const size_t numNonDetectors) {
  BeamlineCacheData data;
  data.detectorIds = std::make_shared<std::vector<detid_t>>();
  data.detectorPositions = std::make_shared<BeamlineCacheData::Vector3ds>();
  data.detectorRotations = std::make_shared<BeamlineCacheData::Quaternions>();
  data.positions = std::make_shared<BeamlineCacheData::Vector3ds>();
  data.rotations = std::make_shared<BeamlineCacheData::Quaternions>();
  for (size_t i = 0; i < numDetectors; ++i) {
    data.detectorIds->emplace_back(static_cast<detid_t>(100 + i));
    data.detectorPositions->emplace_back(Eigen::Vector3d{0.5 * double(i), 1.25, -3.75});
    data.detectorRotations->emplace_back(Eigen::Quaterniond{0.5, 0.5, 0.5, 0.5 + 0.01 * double(i)}.normalized());
  }
  for (size_t i = 0; i < numNonDetectors; ++i) {
    data.positions->emplace_back(Eigen::Vector3d{-1.0, 2.0 * double(i), 7.5});
    data.rotations->emplace_back(Eigen::Quaterniond{1.0, 0.0, 0.25 * double(i), 0.0}.normalized());
  }
  return data;
}

/// Flatten an instrument and return the resulting wrappers, optionally going
/// through a cache file.
std::pair<std::unique_ptr<ComponentInfo>, std::unique_ptr<DetectorInfo>>
flatten(const std::shared_ptr<const Instrument> &instrument, const std::string &cacheFile = std::string()) {
  return InstrumentVisitor::makeWrappers(*instrument, nullptr, cacheFile);
}

/// Write a cache file describing @p instrument, bypassing the detector-count
/// threshold that stops small instruments being cached automatically.
void writeCacheFor(const std::shared_ptr<const Instrument> &instrument, const std::string &cacheFile) {
  InstrumentVisitor visitor(instrument);
  visitor.walkInstrument();
  TS_ASSERT(BeamlineCache::write(cacheFile, visitor.cacheData()));
}

} // namespace

class BeamlineCacheTest : public CxxTest::TestSuite {
public:
  static BeamlineCacheTest *createSuite() { return new BeamlineCacheTest(); }
  static void destroySuite(BeamlineCacheTest *suite) { delete suite; }

  void test_round_trip_preserves_every_value() {
    TemporaryCacheFile file("BeamlineCacheTest_roundtrip.beamline");
    const auto written = makeData(7, 3);
    TS_ASSERT(BeamlineCache::write(file.path(), written));
    TS_ASSERT(file.exists());

    BeamlineCacheData read;
    TS_ASSERT(BeamlineCache::read(file.path(), read));
    TS_ASSERT(read.isComplete());

    TS_ASSERT_EQUALS(*read.detectorIds, *written.detectorIds);
    TS_ASSERT_EQUALS(read.detectorPositions->size(), written.detectorPositions->size());
    for (size_t i = 0; i < written.detectorPositions->size(); ++i) {
      // Bit-for-bit: the point of the cache is to stand in for the values the
      // tree walk would have produced, so any drift is a defect.
      TS_ASSERT_EQUALS((*read.detectorPositions)[i], (*written.detectorPositions)[i]);
      TS_ASSERT_EQUALS((*read.detectorRotations)[i].coeffs(), (*written.detectorRotations)[i].coeffs());
    }
    TS_ASSERT_EQUALS(read.positions->size(), written.positions->size());
    for (size_t i = 0; i < written.positions->size(); ++i) {
      TS_ASSERT_EQUALS((*read.positions)[i], (*written.positions)[i]);
      TS_ASSERT_EQUALS((*read.rotations)[i].coeffs(), (*written.rotations)[i].coeffs());
    }
  }

  void test_round_trip_with_no_non_detector_components() {
    TemporaryCacheFile file("BeamlineCacheTest_detectors_only.beamline");
    TS_ASSERT(BeamlineCache::write(file.path(), makeData(4, 0)));

    BeamlineCacheData read;
    TS_ASSERT(BeamlineCache::read(file.path(), read));
    TS_ASSERT_EQUALS(read.detectorIds->size(), 4);
    TS_ASSERT_EQUALS(read.positions->size(), 0);
  }

  void test_read_of_missing_file_fails_quietly() {
    TemporaryCacheFile file("BeamlineCacheTest_absent.beamline");
    BeamlineCacheData read;
    TS_ASSERT(!BeamlineCache::read(file.path(), read));
    TS_ASSERT(!read.isComplete());
  }

  void test_incomplete_data_is_not_written() {
    TemporaryCacheFile file("BeamlineCacheTest_incomplete.beamline");
    BeamlineCacheData data;
    TS_ASSERT(!BeamlineCache::write(file.path(), data));
    TS_ASSERT(!file.exists());

    // Detector arrays whose lengths disagree are also refused.
    data = makeData(3, 1);
    data.detectorRotations->pop_back();
    TS_ASSERT(!BeamlineCache::write(file.path(), data));
    TS_ASSERT(!file.exists());
  }

  void test_truncated_file_is_rejected_and_deleted() {
    TemporaryCacheFile file("BeamlineCacheTest_truncated.beamline");
    TS_ASSERT(BeamlineCache::write(file.path(), makeData(5, 2)));

    std::filesystem::resize_file(file.path(), std::filesystem::file_size(file.path()) - 16);

    BeamlineCacheData read;
    TS_ASSERT(!BeamlineCache::read(file.path(), read));
    // Deleted, so the next load rewrites it rather than failing the same way.
    TS_ASSERT(!file.exists());
  }

  void test_file_with_foreign_header_is_rejected_and_deleted() {
    TemporaryCacheFile file("BeamlineCacheTest_foreign.beamline");
    {
      std::ofstream stream(file.path(), std::ios::binary);
      const std::vector<char> rubbish(256, '\x7f');
      stream.write(rubbish.data(), static_cast<std::streamsize>(rubbish.size()));
    }

    BeamlineCacheData read;
    TS_ASSERT(!BeamlineCache::read(file.path(), read));
    TS_ASSERT(!file.exists());
  }

  void test_flattening_through_a_cache_matches_flattening_from_the_tree() {
    TemporaryCacheFile file("BeamlineCacheTest_instrument.beamline");
    auto instrument = ComponentCreationHelper::createTestInstrumentRectangular(2 /*banks*/, 4 /*pixels*/);
    writeCacheFor(instrument, file.path());

    const auto direct = flatten(instrument);
    const auto cached = flatten(instrument, file.path());

    assertSameGeometry(*direct.first, *cached.first, *direct.second, *cached.second);
  }

  void test_flattening_through_a_cache_matches_for_a_cylindrical_instrument() {
    TemporaryCacheFile file("BeamlineCacheTest_cylindrical.beamline");
    auto instrument = ComponentCreationHelper::createTestInstrumentCylindrical(3 /*banks*/);
    writeCacheFor(instrument, file.path());

    const auto direct = flatten(instrument);
    const auto cached = flatten(instrument, file.path());

    assertSameGeometry(*direct.first, *cached.first, *direct.second, *cached.second);
  }

  void test_cache_written_for_another_instrument_is_ignored() {
    TemporaryCacheFile file("BeamlineCacheTest_mismatch.beamline");
    auto other = ComponentCreationHelper::createTestInstrumentRectangular(3 /*banks*/, 5 /*pixels*/);
    writeCacheFor(other, file.path());

    // Same file, different instrument: the detector IDs will not line up, so
    // the walk must fall back to deriving the positions from the tree.
    auto instrument = ComponentCreationHelper::createTestInstrumentRectangular(2 /*banks*/, 4 /*pixels*/);
    const auto direct = flatten(instrument);
    const auto cached = flatten(instrument, file.path());

    assertSameGeometry(*direct.first, *cached.first, *direct.second, *cached.second);
    // Dropped, so it is not re-read and rejected on every later load.
    TS_ASSERT(!file.exists());
  }

  void test_a_missing_cache_file_is_written_only_for_large_instruments() {
    TemporaryCacheFile file("BeamlineCacheTest_threshold.beamline");
    auto instrument = ComponentCreationHelper::createTestInstrumentRectangular(2 /*banks*/, 4 /*pixels*/);

    const auto wrappers = flatten(instrument, file.path());
    TS_ASSERT_LESS_THAN(wrappers.second->size(), BeamlineCache::MIN_DETECTORS);
    TS_ASSERT(!file.exists());
  }

private:
  void assertSameGeometry(const ComponentInfo &expectedComponents, const ComponentInfo &actualComponents,
                          const DetectorInfo &expectedDetectors, const DetectorInfo &actualDetectors) {
    TS_ASSERT_EQUALS(expectedComponents.size(), actualComponents.size());
    TS_ASSERT_EQUALS(expectedDetectors.size(), actualDetectors.size());

    for (size_t i = 0; i < expectedComponents.size(); ++i) {
      TS_ASSERT_EQUALS(expectedComponents.position(i), actualComponents.position(i));
      TS_ASSERT_EQUALS(expectedComponents.rotation(i), actualComponents.rotation(i));
    }
    for (size_t i = 0; i < expectedDetectors.size(); ++i) {
      TS_ASSERT_EQUALS(expectedDetectors.position(i), actualDetectors.position(i));
      TS_ASSERT_EQUALS(expectedDetectors.rotation(i), actualDetectors.rotation(i));
      TS_ASSERT_EQUALS(expectedDetectors.isMonitor(i), actualDetectors.isMonitor(i));
    }
    TS_ASSERT_EQUALS(expectedComponents.source(), actualComponents.source());
    TS_ASSERT_EQUALS(expectedComponents.sample(), actualComponents.sample());
  }
};
