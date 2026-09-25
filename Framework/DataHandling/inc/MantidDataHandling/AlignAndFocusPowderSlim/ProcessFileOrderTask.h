// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +

#pragma once

#include "MantidAPI/Progress.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/BankCalibration.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/DirectEventReader.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/NexusLoader.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/ProcessBankTaskBase.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/SpectraProcessingData.h"

#include <H5Cpp.h>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace Mantid::DataHandling::AlignAndFocusPowderSlim {

/** Reads the events of all banks in one pass through the file, in file order, and histograms them.
 *
 * The SNS data acquisition writes the event columns of all banks interleaved, in runs of a few MB. Reading bank by
 * bank, or column by column, asks the file system for blocks that are mostly other columns' data, and on network
 * storage that reads large blocks this multiplies the bytes fetched. Here the events of every bank are cut into
 * batches, the batches are ordered by where they sit in the file and grouped into waves, and each wave is read as
 * block-aligned pieces in file order on a pool of threads, so each block is fetched once. The next wave is read while
 * the current one is histogrammed.
 *
 * Histogramming is the same as ProcessBankTask's; counts are sums, so the order in which batches arrive does not
 * change the result. Only used without splitters, and only when every bank's columns can be read directly.
 */
class ProcessFileOrderTask : public ProcessBankTaskBase {
public:
  /// Plans the batches and waves; call run() to read and histogram them.
  ProcessFileOrderTask(std::vector<std::string> &bankEntryNames, H5::H5File &h5file,
                       std::shared_ptr<NexusLoader> loader, std::shared_ptr<const DirectEventReader> reader,
                       SpectraProcessingData &processingData, const BankCalibrationFactory &calibFactory,
                       const size_t events_per_chunk, const size_t grainsize_event);

  /// True when every bank's detector-id and time-of-flight columns can be read directly.
  static bool canProcess(const DirectEventReader &reader, const std::vector<std::string> &bankEntryNames);

  /// Number of waves, one progress step each.
  size_t numWaves() const { return m_waves.size(); }

  void run(API::Progress &progress);

private:
  struct Batch {
    size_t bank;
    std::vector<size_t> offsets;
    std::vector<size_t> slabsizes;
    size_t numEvents{0};
    uint64_t filePosition{0};
    /// where the batch's events land, in the buffer of the wave it belongs to
    std::span<uint32_t> detid;
    std::span<float> tof;
  };
  /// Room for one wave's events. Two are allocated once and alternate, so that one wave is read while the previous
  /// one is histogrammed; allocating and zero-filling buffers for every wave took as long as the histogramming.
  struct WaveBuffer {
    std::unique_ptr<uint32_t[]> detid;
    std::unique_ptr<float[]> tof;
  };
  struct Bank {
    std::string detidPath;
    std::string tofPath;
    std::vector<BankCalibration> calibrations;
  };
  /// A wave is batches [first, last) of m_batches.
  using Wave = std::pair<size_t, size_t>;

  std::vector<ParallelFileReader::Task> prepareWave(const Wave &wave, WaveBuffer &buffer);
  void histogramWave(const Wave &wave);

  std::shared_ptr<const DirectEventReader> m_reader;
  SpectraProcessingData &m_processingData;
  const size_t m_grainsize_event;
  std::vector<Bank> m_banks;
  std::vector<Batch> m_batches;
  std::vector<Wave> m_waves;
  /// the most events in any wave
  size_t m_waveCapacity{0};
};

} // namespace Mantid::DataHandling::AlignAndFocusPowderSlim
