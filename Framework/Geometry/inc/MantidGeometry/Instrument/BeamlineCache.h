// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

#include "MantidGeometry/DllConfig.h"
#include "MantidGeometry/IDTypes.h"

#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace Mantid {
namespace Geometry {

/** The bulk numeric arrays of the flattened instrument, in the struct-of-arrays
 * layout that Beamline::ComponentInfo and Beamline::DetectorInfo consume.
 *
 * Only the arrays that are both expensive to derive from the component tree and
 * large enough to be worth storing are held here. Absolute positions and
 * rotations qualify: Component::getPos() re-walks the parent chain and calls
 * getRotation() at each level, so deriving them is quadratic in the depth of
 * the tree and is repeated once per detector. Everything else the tree walk
 * produces - names, scale factors, component types, the assembly ranges and the
 * child lists - is either a direct member read or is sized by the number of
 * non-detector components, which is small, so it is recomputed rather than
 * stored.
 *
 * The arrays are held by shared_ptr because both the producer
 * (InstrumentVisitor) and the consumer own them through shared_ptr already, so
 * a cache round trip copies no geometry.
 */
struct MANTID_GEOMETRY_DLL BeamlineCacheData {
  using Vector3ds = std::vector<Eigen::Vector3d>;
  using Quaternions = std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>>;

  /// Detector IDs in visitation order. Stored so that a cache file can be
  /// matched against the instrument it is about to be applied to.
  std::shared_ptr<std::vector<detid_t>> detectorIds;
  /// Absolute positions of the detectors, indexed by detector index
  std::shared_ptr<Vector3ds> detectorPositions;
  /// Absolute rotations of the detectors, indexed by detector index
  std::shared_ptr<Quaternions> detectorRotations;
  /// Absolute positions of the non-detector components, indexed by component
  /// index minus the number of detectors
  std::shared_ptr<Vector3ds> positions;
  /// Absolute rotations of the non-detector components, indexed as above
  std::shared_ptr<Quaternions> rotations;

  bool isComplete() const;
};

/** BeamlineCache : reads and writes the flattened instrument's position and
 * rotation arrays to a file next to the `.vtp` geometry cache.
 *
 * The `.vtp` cache holds tessellated CSG shapes only, so the flattening of the
 * component tree into struct-of-arrays form is repeated on every load. This
 * cache covers the expensive part of that flattening.
 *
 * The file is a fixed-width, 16-byte-aligned, struct-of-arrays dump with a
 * 64-byte header. There are no variable-length records, so it can be consumed
 * by a memory map rather than a read once Beamline::ComponentInfo is able to
 * reference memory it does not own. It is keyed the same way as the `.vtp`
 * file, on the instrument name plus a checksum of the IDF text, and it is
 * machine-local: a file written by a build with a different floating-point
 * layout or byte order is rejected and regenerated rather than converted.
 */
class MANTID_GEOMETRY_DLL BeamlineCache {
public:
  /// Bumped whenever the on-disk layout changes. Files written by an older
  /// version are rejected and rewritten.
  static constexpr uint32_t FORMAT_VERSION = 1;

  /// Instruments with fewer detectors than this are not worth a cache file: the
  /// tree walk they avoid is shorter than the time spent opening the file, and
  /// the 506-file instrument corpus would otherwise fill the cache directory
  /// with files that save nothing.
  static constexpr size_t MIN_DETECTORS = 100000;

  /// Suffix used for the cache file, appended to the mangled instrument name
  static const std::string FILE_EXTENSION;

  /// Read @p filename into @p data. Returns false, and removes the file, if it
  /// is absent, truncated, or was written by an incompatible build.
  static bool read(const std::string &filename, BeamlineCacheData &data);

  /// Write @p data to @p filename, via a temporary file so that an interrupted
  /// write cannot leave a partial cache behind. Returns false on any failure;
  /// failing to write a cache is never an error for the caller.
  static bool write(const std::string &filename, const BeamlineCacheData &data);

  /// Delete a cache file that turned out not to describe the instrument it was
  /// read for, so that it is not read again on every subsequent load.
  static void discard(const std::string &filename);
};

} // namespace Geometry
} // namespace Mantid
