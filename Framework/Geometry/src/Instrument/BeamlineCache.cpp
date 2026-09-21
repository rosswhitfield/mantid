// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidGeometry/Instrument/BeamlineCache.h"
#include "MantidKernel/Logger.h"

#include <Poco/Process.h>

#include <array>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>

namespace Mantid::Geometry {
namespace {
Kernel::Logger g_log("BeamlineCache");

constexpr std::array<char, 8> MAGIC{{'M', 'T', 'D', 'B', 'E', 'A', 'M', '\0'}};
constexpr size_t HEADER_SIZE = 64;
constexpr size_t ALIGNMENT = 16;

// The arrays are written as raw memory, so the reader has to be sure that the
// writer laid them out the same way. Both are fixed by Eigen, but assert it
// rather than assume it.
static_assert(sizeof(Eigen::Vector3d) == 3 * sizeof(double), "Eigen::Vector3d must be three packed doubles");
static_assert(sizeof(Eigen::Quaterniond) == 4 * sizeof(double), "Eigen::Quaterniond must be four packed doubles");

/** A value that differs between builds whose raw arrays are not interchangeable.
 * A cache file is only ever read back by the machine that wrote it, so rather
 * than defining a portable encoding we record enough of the local layout to
 * recognise a file that must be discarded. */
uint64_t layoutMarker() {
  // 0x0102030405060708 read back byte by byte distinguishes the byte order, and
  // the sizes cover the rest of what the raw dump depends on.
  constexpr uint64_t probe = 0x0102030405060708ULL;
  std::array<unsigned char, sizeof(probe)> bytes{};
  std::memcpy(bytes.data(), &probe, sizeof(probe));
  uint64_t marker = 0;
  for (const auto byte : bytes)
    marker = (marker << 8) | byte;
  return marker ^ (static_cast<uint64_t>(sizeof(double)) << 32) ^ (static_cast<uint64_t>(sizeof(size_t)) << 40) ^
         (static_cast<uint64_t>(sizeof(detid_t)) << 48);
}

/// Round @p bytes up to the next multiple of ALIGNMENT, so that every array in
/// the file starts on an address a memory map could hand straight to Eigen.
constexpr size_t padded(const size_t bytes) { return ((bytes + ALIGNMENT - 1) / ALIGNMENT) * ALIGNMENT; }

struct Header {
  std::array<char, 8> magic{};
  uint32_t version{0};
  uint32_t reserved{0};
  uint64_t marker{0};
  uint64_t numDetectors{0};
  uint64_t numNonDetectors{0};
  uint64_t payloadBytes{0};
};
static_assert(sizeof(Header) <= HEADER_SIZE, "Header must fit in the reserved space");

/// Byte length of the payload that follows the header, with each array padded
/// up to the alignment boundary.
size_t payloadSize(const size_t numDetectors, const size_t numNonDetectors) {
  return padded(numDetectors * sizeof(detid_t)) + padded(numDetectors * sizeof(Eigen::Vector3d)) +
         padded(numDetectors * sizeof(Eigen::Quaterniond)) + padded(numNonDetectors * sizeof(Eigen::Vector3d)) +
         padded(numNonDetectors * sizeof(Eigen::Quaterniond));
}

/// Write @p vec as raw bytes followed by enough zeros to reach the next
/// alignment boundary. Templated on the vector rather than its element so that
/// the Eigen-aligned allocator is handled by the same overload.
template <typename Vector> void writeArray(std::ostream &os, const Vector &vec) {
  const size_t bytes = vec.size() * sizeof(typename Vector::value_type);
  if (bytes > 0)
    os.write(reinterpret_cast<const char *>(vec.data()), static_cast<std::streamsize>(bytes));
  static const std::array<char, ALIGNMENT> zeros{};
  os.write(zeros.data(), static_cast<std::streamsize>(padded(bytes) - bytes));
}

/// Read @p count elements into a freshly sized vector, then skip the padding.
template <typename Vector> bool readArray(std::istream &is, Vector &vec, const size_t count) {
  vec.resize(count);
  const size_t bytes = count * sizeof(typename Vector::value_type);
  if (bytes > 0 && !is.read(reinterpret_cast<char *>(vec.data()), static_cast<std::streamsize>(bytes)))
    return false;
  return static_cast<bool>(is.seekg(static_cast<std::streamoff>(padded(bytes) - bytes), std::ios::cur));
}

} // namespace

const std::string BeamlineCache::FILE_EXTENSION = ".beamline";

void BeamlineCache::discard(const std::string &filename) {
  std::error_code ec;
  std::filesystem::remove(filename, ec);
}

bool BeamlineCacheData::isComplete() const {
  return detectorIds && detectorPositions && detectorRotations && positions && rotations &&
         detectorIds->size() == detectorPositions->size() && detectorIds->size() == detectorRotations->size() &&
         positions->size() == rotations->size();
}

bool BeamlineCache::read(const std::string &filename, BeamlineCacheData &data) {
  std::error_code ec;
  if (!std::filesystem::exists(filename, ec))
    return false;

  std::ifstream file(filename, std::ios::binary);
  if (!file)
    return false;

  Header header;
  if (!file.read(reinterpret_cast<char *>(&header), sizeof(header))) {
    g_log.information() << "Beamline cache " << filename << " is too short to hold a header; ignoring it.\n";
    BeamlineCache::discard(filename);
    return false;
  }
  if (header.magic != MAGIC || header.version != FORMAT_VERSION || header.marker != layoutMarker()) {
    g_log.information() << "Beamline cache " << filename << " was written by an incompatible build; ignoring it.\n";
    BeamlineCache::discard(filename);
    return false;
  }
  if (header.payloadBytes != payloadSize(header.numDetectors, header.numNonDetectors)) {
    g_log.information() << "Beamline cache " << filename << " has an inconsistent header; ignoring it.\n";
    BeamlineCache::discard(filename);
    return false;
  }
  const auto fileSize = std::filesystem::file_size(filename, ec);
  if (ec || fileSize != HEADER_SIZE + header.payloadBytes) {
    g_log.information() << "Beamline cache " << filename << " is truncated; ignoring it.\n";
    BeamlineCache::discard(filename);
    return false;
  }
  if (!file.seekg(HEADER_SIZE, std::ios::beg)) {
    BeamlineCache::discard(filename);
    return false;
  }

  auto detectorIds = std::make_shared<std::vector<detid_t>>();
  auto detectorPositions = std::make_shared<BeamlineCacheData::Vector3ds>();
  auto detectorRotations = std::make_shared<BeamlineCacheData::Quaternions>();
  auto positions = std::make_shared<BeamlineCacheData::Vector3ds>();
  auto rotations = std::make_shared<BeamlineCacheData::Quaternions>();

  const auto numDetectors = static_cast<size_t>(header.numDetectors);
  const auto numNonDetectors = static_cast<size_t>(header.numNonDetectors);
  if (!readArray(file, *detectorIds, numDetectors) || !readArray(file, *detectorPositions, numDetectors) ||
      !readArray(file, *detectorRotations, numDetectors) || !readArray(file, *positions, numNonDetectors) ||
      !readArray(file, *rotations, numNonDetectors)) {
    g_log.information() << "Beamline cache " << filename << " could not be read; ignoring it.\n";
    BeamlineCache::discard(filename);
    return false;
  }

  data.detectorIds = std::move(detectorIds);
  data.detectorPositions = std::move(detectorPositions);
  data.detectorRotations = std::move(detectorRotations);
  data.positions = std::move(positions);
  data.rotations = std::move(rotations);
  return true;
}

bool BeamlineCache::write(const std::string &filename, const BeamlineCacheData &data) {
  if (!data.isComplete())
    return false;

  // Write to a sibling temporary and rename, so that a crash cannot leave a
  // half-written file that looks valid. The name carries the process id
  // because several processes may cache the same instrument at once.
  const auto temporary = filename + ".tmp" + std::to_string(Poco::Process::id());
  {
    std::ofstream file(temporary, std::ios::binary | std::ios::trunc);
    if (!file) {
      g_log.information() << "Cannot write the beamline cache to " << temporary << "; it will not be cached.\n";
      return false;
    }

    Header header;
    header.magic = MAGIC;
    header.version = FORMAT_VERSION;
    header.marker = layoutMarker();
    header.numDetectors = data.detectorIds->size();
    header.numNonDetectors = data.positions->size();
    header.payloadBytes = payloadSize(data.detectorIds->size(), data.positions->size());

    std::array<char, HEADER_SIZE> headerBlock{};
    std::memcpy(headerBlock.data(), &header, sizeof(header));
    file.write(headerBlock.data(), HEADER_SIZE);

    writeArray(file, *data.detectorIds);
    writeArray(file, *data.detectorPositions);
    writeArray(file, *data.detectorRotations);
    writeArray(file, *data.positions);
    writeArray(file, *data.rotations);

    if (!file) {
      g_log.information() << "Failed while writing the beamline cache " << temporary << ".\n";
      file.close();
      BeamlineCache::discard(temporary);
      return false;
    }
  }

  std::error_code ec;
  std::filesystem::rename(temporary, filename, ec);
  if (ec) {
    g_log.information() << "Cannot move the beamline cache into place at " << filename << ": " << ec.message() << "\n";
    BeamlineCache::discard(temporary);
    return false;
  }
  return true;
}

} // namespace Mantid::Geometry
