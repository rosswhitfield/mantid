// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidDataHandling/LoadEventNexusMD.h"
#include "MantidDataHandling/AlignAndFocusPowderSlim/NexusLoader.h"
#include "MantidDataHandling/LoadEventNexus.h"
#include "MantidDataHandling/LoadEventNexusMD/ProcessEventsMDTask.h"

#include "MantidAPI/ExperimentInfo.h"
#include "MantidAPI/FileProperty.h"
#include "MantidAPI/IMDEventWorkspace.h"
#include "MantidAPI/Progress.h"
#include "MantidAPI/Run.h"
#include "MantidAPI/Sample.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidDataObjects/MDEventFactory.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidGeometry/Instrument.h"
#include "MantidGeometry/Instrument/DetectorInfo.h"
#include "MantidGeometry/Instrument/Goniometer.h"
#include "MantidGeometry/MDGeometry/MDFrameFactory.h"
#include "MantidGeometry/MDGeometry/MDHistoDimension.h"
#include "MantidGeometry/MDGeometry/QLab.h"
#include "MantidGeometry/MDGeometry/QSample.h"
#include "MantidKernel/ArrayProperty.h"
#include "MantidKernel/BoundedValidator.h"
#include "MantidKernel/ConfigService.h"
#include "MantidKernel/ListValidator.h"
#include "MantidKernel/ThreadPool.h"
#include "MantidKernel/ThreadScheduler.h"
#include "MantidKernel/TimeSeriesProperty.h"
#include "MantidKernel/Unit.h"
#include "MantidNexus/H5Util.h"

#include <algorithm>
#include <limits>
#include <regex>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_reduce.h>

namespace Mantid::DataHandling::LoadEventNexusMD {

using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::Kernel;
using Mantid::DataHandling::AlignAndFocusPowderSlim::NexusLoader;
namespace NxsFieldNames = Mantid::DataHandling::AlignAndFocusPowderSlim::NxsFieldNames;

namespace {
/// NeXus entry that everything else hangs off
const std::string ENTRY_TOP_LEVEL("entry");
const std::string MICROSEC("microseconds");
/// grown/shrunk from the calculated extents so events exactly on a boundary are not rejected
constexpr double EXTENTS_PADDING{1.e-5};
/**
 * Calculated extents are clamped to this, in inverse Angstrom. No neutron source reaches it - it corresponds to a
 * wavelength of 0.06 Angstrom - but a stray event with a time-of-flight of a microsecond or two would otherwise blow
 * the box up and waste the whole recursion depth on empty space.
 */
constexpr double MAX_CALCULATED_Q{200.};
Logger g_log("LoadEventNexusMD");

/**
 * tbb::parallel_reduce body finding the smallest strictly positive corrected time-of-flight, and the largest, in a
 * block of raw values. Non-positive times cannot produce a finite Q, so they must not be allowed to set the range.
 */
struct TofRangeTask {
  TofRangeTask(const std::vector<float> *tofs, const double scale, const double offset)
      : m_tofs(tofs), m_scale(scale), m_offset(offset) {}
  TofRangeTask(const TofRangeTask &other, tbb::split)
      : m_tofs(other.m_tofs), m_scale(other.m_scale), m_offset(other.m_offset) {}

  void operator()(const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      const double tof = static_cast<double>((*m_tofs)[i]) * m_scale + m_offset;
      if (tof <= 0.)
        continue;
      minPositive = std::min(minPositive, tof);
      maxValue = std::max(maxValue, tof);
    }
  }

  void join(const TofRangeTask &other) {
    minPositive = std::min(minPositive, other.minPositive);
    maxValue = std::max(maxValue, other.maxValue);
  }

  double minPositive{std::numeric_limits<double>::max()};
  double maxValue{std::numeric_limits<double>::lowest()};

private:
  const std::vector<float> *m_tofs;
  double m_scale;
  double m_offset;
};
} // namespace

DECLARE_ALGORITHM(LoadEventNexusMD)

const std::string LoadEventNexusMD::name() const { return "LoadEventNexusMD"; }

int LoadEventNexusMD::version() const { return 1; }

const std::string LoadEventNexusMD::category() const { return "DataHandling\\Nexus;MDAlgorithms\\Creation"; }

const std::string LoadEventNexusMD::summary() const {
  return "Load single crystal diffraction event data from a NeXus file directly into an MDEventWorkspace in Q.";
}

const std::vector<std::string> LoadEventNexusMD::seeAlso() const {
  return {"LoadEventNexus", "ConvertToMD", "ConvertToDiffractionMDWorkspace", "AlignAndFocusPowderSlim"};
}

void LoadEventNexusMD::init() {
  const std::vector<std::string> exts{".nxs.h5", ".nxs", "_event.nxs"};
  declareProperty(std::make_unique<FileProperty>(PropertyNames::FILENAME, "", FileProperty::Load, exts),
                  "The name of the Event NeXus file to read, including its full or relative path.");

  declareProperty(
      std::make_unique<WorkspaceProperty<IMDEventWorkspace>>(PropertyNames::OUTPUT_WKSP, "", Direction::Output),
      "The name of the output MDEventWorkspace.");

  declareProperty(PropertyNames::Q_FRAME, QFrames::SAMPLE,
                  std::make_shared<StringListValidator>(std::vector<std::string>{QFrames::SAMPLE, QFrames::LAB}),
                  "Frame of the output dimensions.\n"
                  "  " +
                      QFrames::SAMPLE +
                      ": wave-vector change in the frame of the sample, i.e. with the "
                      "goniometer rotation taken out.\n"
                      "  " +
                      QFrames::LAB + ": wave-vector change in the laboratory frame.");

  declareProperty(PropertyNames::EVENT_TYPE, EventTypes::FULL,
                  std::make_shared<StringListValidator>(std::vector<std::string>{EventTypes::FULL, EventTypes::LEAN}),
                  "Type of event stored in the output workspace.\n"
                  "  " +
                      EventTypes::FULL +
                      ": keeps the detector id and goniometer index of every event, which "
                      "peak integration and normalisation need.\n"
                      "  " +
                      EventTypes::LEAN + ": signal, error and coordinates only, using about half the memory.");

  declareProperty(std::make_unique<ArrayProperty<double>>(PropertyNames::MIN_VALUES),
                  "Three comma separated values, the smallest Q in each dimension. Events below these values are "
                  "not added to the workspace. Leave empty to calculate the enclosing box from the instrument "
                  "geometry and the time-of-flight range of the file.");
  declareProperty(std::make_unique<ArrayProperty<double>>(PropertyNames::MAX_VALUES),
                  "Three comma separated values, the largest Q in each dimension. Events at or above these values "
                  "are not added to the workspace. Must be given together with " +
                      PropertyNames::MIN_VALUES + ".");

  // Box controller properties, matching ConvertToMD. SplitInto matters more than it looks: every event walks the box
  // tree from the root, so the shallower tree that a larger SplitInto gives is substantially faster to fill. On a 128
  // million event TOPAZ run, splitting into 2 rather than 5 costs about 70% more time.
  this->initBoxControllerProps("5" /*SplitInto*/, 1000 /*SplitThreshold*/, 20 /*MaxRecursionDepth*/);

  auto mustBePositive = std::make_shared<BoundedValidator<int>>();
  mustBePositive->setLower(1);
  declareProperty(PropertyNames::READ_SIZE_FROM_DISK, 10000000, mustBePositive,
                  "Number of events to read from the file at a time.");
  declareProperty(PropertyNames::EVENTS_PER_THREAD, 1000, mustBePositive,
                  "Number of events to process in a single thread. Higher values reduce threading overhead.");
  const std::string chunkingGroup("Chunking");
  setPropertyGroup(PropertyNames::READ_SIZE_FROM_DISK, chunkingGroup);
  setPropertyGroup(PropertyNames::EVENTS_PER_THREAD, chunkingGroup);
}

std::map<std::string, std::string> LoadEventNexusMD::validateInputs() {
  std::map<std::string, std::string> errors;

  const int readSize = getProperty(PropertyNames::READ_SIZE_FROM_DISK);
  const int perThread = getProperty(PropertyNames::EVENTS_PER_THREAD);
  if (readSize < perThread)
    errors[PropertyNames::READ_SIZE_FROM_DISK] =
        PropertyNames::READ_SIZE_FROM_DISK + " must not be smaller than " + PropertyNames::EVENTS_PER_THREAD;

  const std::vector<double> minVals = getProperty(PropertyNames::MIN_VALUES);
  const std::vector<double> maxVals = getProperty(PropertyNames::MAX_VALUES);
  if (minVals.empty() != maxVals.empty()) {
    const std::string msg(PropertyNames::MIN_VALUES + " and " + PropertyNames::MAX_VALUES +
                          " must either both be set or both be left empty");
    errors[PropertyNames::MIN_VALUES] = msg;
    errors[PropertyNames::MAX_VALUES] = msg;
  } else if (!minVals.empty()) {
    if (minVals.size() != 3 || maxVals.size() != 3) {
      const std::string msg("must have 3 values, one for each dimension");
      errors[PropertyNames::MIN_VALUES] = msg;
      errors[PropertyNames::MAX_VALUES] = msg;
    } else {
      for (size_t i = 0; i < 3; ++i) {
        if (!(minVals[i] < maxVals[i])) {
          const std::string msg("max is not greater than min at index " + std::to_string(i));
          errors[PropertyNames::MIN_VALUES] = msg;
          errors[PropertyNames::MAX_VALUES] = msg;
          break;
        }
      }
    }
  }

  return errors;
}

/**
 * Find every NXevent_data entry directly below /entry. Note that the entry name need not correspond to an instrument
 * component: the detectors are found from the event_id of each event instead.
 */
std::vector<std::string>
LoadEventNexusMD::determineBanksToLoad(const Mantid::Nexus::NexusDescriptor &descriptor) const {
  const std::set<std::string> classEntries = descriptor.allAddressesOfType("NXevent_data");
  if (classEntries.empty())
    throw std::runtime_error("No NXevent_data entries found in file");

  std::vector<std::string> bankEntryNames;
  const std::regex classRegex("(/" + ENTRY_TOP_LEVEL + "/)([^/]*)");
  std::smatch groups;
  for (const std::string &classEntry : classEntries) {
    if (!std::regex_match(classEntry, groups, classRegex))
      continue;
    if (classEntry.ends_with("bank_error_events") || classEntry.ends_with("bank_unmapped_events"))
      continue;
    bankEntryNames.push_back(groups[2].str());
  }
  return bankEntryNames;
}

/**
 * Single crystal instruments carry an additive time-of-flight correction in their parameter file. LoadEventNexus
 * applies it to every event and records it as the T0 run property, so this algorithm has to do the same or the
 * wavelengths, and therefore every Q, come out slightly wrong.
 */
double LoadEventNexusMD::getT0(const MatrixWorkspace_sptr &wksp) const {
  const auto instrument = wksp->getInstrument();
  if (!instrument->hasParameter("T0"))
    return 0.;
  const std::vector<double> t0 = instrument->getNumberParameter("T0", true);
  return t0.empty() ? 0. : t0.front();
}

std::pair<double, double> LoadEventNexusMD::findTimeOfFlightRange(H5::H5File &h5file,
                                                                  const std::vector<std::string> &bankEntryNames,
                                                                  const double t0) const {
  const auto events_per_chunk = static_cast<size_t>(int(getProperty(PropertyNames::READ_SIZE_FROM_DISK)));
  const NexusLoader loader(false, {});

  double tof_min{std::numeric_limits<double>::max()};
  double tof_max{std::numeric_limits<double>::lowest()};

  auto entry = h5file.openGroup(ENTRY_TOP_LEVEL);
  auto tof_values = std::make_unique<std::vector<float>>();
  for (const auto &bankEntryName : bankEntryNames) {
    auto event_group = entry.openGroup(bankEntryName);
    auto tof_SDS = event_group.openDataSet(NxsFieldNames::TIME_OF_FLIGHT);
    const auto total_events = static_cast<size_t>(tof_SDS.getSpace().getSelectNpoints());
    if (total_events == 0)
      continue;

    std::string tof_unit;
    Nexus::H5Util::readStringAttribute(tof_SDS, "units", tof_unit);
    const double time_conversion = Units::timeConversionValue(tof_unit, MICROSEC);

    for (size_t offset = 0; offset < total_events; offset += events_per_chunk) {
      const std::vector<size_t> offsets{offset};
      const std::vector<size_t> slabsizes{std::min(events_per_chunk, total_events - offset)};
      loader.loadData(tof_SDS, tof_values, offsets, slabsizes);
      if (tof_values->empty())
        continue;
      TofRangeTask task(tof_values.get(), time_conversion, t0);
      tbb::parallel_reduce(tbb::blocked_range<size_t>(0, tof_values->size()), task);
      tof_min = std::min(tof_min, task.minPositive);
      tof_max = std::max(tof_max, task.maxValue);
    }
  }

  return {tof_min, tof_max};
}

/**
 * Every detector sees Q along a fixed ray, so once the time-of-flight range is known the enclosing box follows
 * exactly: sweep each direction between the two wavenumbers it can produce and take the componentwise extremes.
 */
void LoadEventNexusMD::calculateExtents(H5::H5File &h5file, const std::vector<std::string> &bankEntryNames,
                                        const QConversionTable &table, const double t0, std::vector<double> &minVals,
                                        std::vector<double> &maxVals) const {
  const auto [tof_min, tof_max] = this->findTimeOfFlightRange(h5file, bankEntryNames, t0);
  if (!(tof_min > 0.) || !(tof_max >= tof_min))
    throw std::runtime_error("Failed to find a usable time-of-flight range in the file. Supply " +
                             PropertyNames::MIN_VALUES + " and " + PropertyNames::MAX_VALUES + " instead.");
  g_log.information() << "Time-of-flight range " << tof_min << " to " << tof_max << " microseconds\n";

  minVals.assign(3, std::numeric_limits<double>::max());
  maxVals.assign(3, std::numeric_limits<double>::lowest());

  for (const auto &coeff : table.coefficients()) {
    if (coeff.kfactor <= 0.)
      continue;
    // the largest wavenumber comes from the shortest time-of-flight
    const double k_hi = coeff.kfactor / tof_min;
    const double k_lo = coeff.kfactor / tof_max;
    const double direction[3] = {coeff.qx, coeff.qy, coeff.qz};
    for (size_t d = 0; d < 3; ++d) {
      const double q_lo = direction[d] * k_lo;
      const double q_hi = direction[d] * k_hi;
      minVals[d] = std::min(minVals[d], std::min(q_lo, q_hi));
      maxVals[d] = std::max(maxVals[d], std::max(q_lo, q_hi));
    }
  }

  bool clamped{false};
  for (size_t d = 0; d < 3; ++d) {
    if (!(minVals[d] < maxVals[d]))
      throw std::runtime_error("Failed to calculate the extents of the output workspace. Supply " +
                               PropertyNames::MIN_VALUES + " and " + PropertyNames::MAX_VALUES + " instead.");
    minVals[d] -= EXTENTS_PADDING;
    maxVals[d] += EXTENTS_PADDING;
    clamped |= (minVals[d] < -MAX_CALCULATED_Q) || (maxVals[d] > MAX_CALCULATED_Q);
    minVals[d] = std::max(minVals[d], -MAX_CALCULATED_Q);
    maxVals[d] = std::min(maxVals[d], MAX_CALCULATED_Q);
  }
  if (clamped)
    g_log.warning() << "Calculated extents were clamped to +/-" << MAX_CALCULATED_Q
                    << " inverse Angstrom, which suggests the file holds events with a very short time-of-flight. "
                       "Events outside the clamped box are not added to the workspace; supply "
                    << PropertyNames::MIN_VALUES << " and " << PropertyNames::MAX_VALUES << " to choose the box.\n";
}

template <typename MDEW>
void LoadEventNexusMD::loadEvents(const std::shared_ptr<MDEW> &ws, H5::H5File &h5file,
                                  const std::vector<std::string> &bankEntryNames, const QConversionTable &table,
                                  const double t0, const Mantid::coord_t *extentsMin,
                                  const Mantid::coord_t *extentsMax) {
  const auto events_per_chunk = static_cast<size_t>(int(getProperty(PropertyNames::READ_SIZE_FROM_DISK)));
  const auto grainsize = static_cast<size_t>(int(getProperty(PropertyNames::EVENTS_PER_THREAD)));
  const NexusLoader loader(false, {});

  auto bc = ws->getBoxController();
  // The scheduler is owned by the pool. Splitting happens between chunks, never while events are being added.
  auto *ts = new ThreadSchedulerFIFO();
  ThreadPool tp(ts);

  size_t nEventsInWS{0};
  size_t eventsAddedSinceSplit{0};
  size_t lastNumBoxes = bc->getTotalNumMDBoxes();

  auto entry = h5file.openGroup(ENTRY_TOP_LEVEL);
  Progress progress(this, .2, .95, bankEntryNames.size());

  // declared once so the memory can be reused between chunks
  auto event_detid = std::make_unique<std::vector<uint32_t>>();
  auto event_tof = std::make_unique<std::vector<float>>();

  for (const auto &bankEntryName : bankEntryNames) {
    auto event_group = entry.openGroup(bankEntryName);
    auto tof_SDS = event_group.openDataSet(NxsFieldNames::TIME_OF_FLIGHT);
    const auto total_events = static_cast<size_t>(tof_SDS.getSpace().getSelectNpoints());
    if (total_events == 0) {
      g_log.debug() << bankEntryName << " empty\n";
      progress.report();
      continue;
    }
    auto detID_SDS = event_group.openDataSet(NxsFieldNames::DETID);

    std::string tof_unit;
    Nexus::H5Util::readStringAttribute(tof_SDS, "units", tof_unit);
    const double tofScale = Units::timeConversionValue(tof_unit, MICROSEC);

    for (size_t offset = 0; offset < total_events; offset += events_per_chunk) {
      const std::vector<size_t> offsets{offset};
      const std::vector<size_t> slabsizes{std::min(events_per_chunk, total_events - offset)};
      tbb::parallel_invoke([&] { loader.loadData(detID_SDS, event_detid, offsets, slabsizes); },
                           [&] { loader.loadData(tof_SDS, event_tof, offsets, slabsizes); });

      std::atomic<size_t> added{0};
      const ProcessEventsMDTask<MDEW> task(event_detid.get(), event_tof.get(), &table, ws, extentsMin, extentsMax,
                                           tofScale, t0, &added, &m_eventsWithNoDetector);
      tbb::parallel_for(tbb::blocked_range<size_t>(0, event_tof->size(), grainsize), task);

      nEventsInWS += added;
      eventsAddedSinceSplit += added;
      if (bc->shouldSplitBoxes(nEventsInWS, eventsAddedSinceSplit, lastNumBoxes)) {
        ws->splitAllIfNeeded(ts);
        tp.joinAll();
        lastNumBoxes = bc->getTotalNumMDBoxes();
        eventsAddedSinceSplit = 0;
      }
    }
    progress.report();
  }

  ws->splitAllIfNeeded(ts);
  tp.joinAll();
  ws->refreshCache();
  g_log.information() << "Added " << nEventsInWS << " events to the output workspace\n";
}

void LoadEventNexusMD::exec() {
  const std::string filename = getPropertyValue(PropertyNames::FILENAME);
  const bool qSample = (getPropertyValue(PropertyNames::Q_FRAME) == QFrames::SAMPLE);
  const std::string eventType = getPropertyValue(PropertyNames::EVENT_TYPE);
  m_eventsWithNoDetector = 0;

  this->progress(.0, "Reading file structure");
  const Mantid::Nexus::NexusDescriptor descriptor(filename);
  const auto bankEntryNames = this->determineBanksToLoad(descriptor);
  g_log.debug() << "Total banks to read: " << bankEntryNames.size() << "\n";

  // A bare Workspace2D is enough to carry the logs and the instrument. It is never initialized with spectra - only
  // DetectorInfo is needed - and ends up as the output workspace's ExperimentInfo.
  this->progress(.02, "Loading logs and instrument");
  MatrixWorkspace_sptr wksp = std::make_shared<Workspace2D>();
  try {
    LoadEventNexus::loadEntryMetadata(filename, wksp, ENTRY_TOP_LEVEL);
  } catch (std::exception &e) {
    g_log.warning() << "Error while loading meta data: " << e.what() << '\n';
  }
  // runLoadNexusLogs also creates a universal goniometer from the omega/chi/phi logs
  auto periodLog = std::make_unique<const TimeSeriesProperty<int>>("period_log"); // not used
  int nPeriods{1};
  LoadEventNexus::runLoadNexusLogs<MatrixWorkspace_sptr>(filename, wksp, *this, false, nPeriods, periodLog);
  LoadEventNexus::loadInstrument<MatrixWorkspace_sptr>(filename, wksp, ENTRY_TOP_LEVEL, this, &descriptor);

  this->progress(.05, "Calculating detector geometry");
  // Q_lab needs no rotation; Q_sample takes the goniometer rotation back out
  DblMatrix rotation(3, 3, true);
  if (qSample) {
    if (!wksp->run().getGoniometer().isDefined())
      g_log.warning() << "No goniometer is defined for this run, so " << QFrames::SAMPLE << " is the same as "
                      << QFrames::LAB << "\n";
    rotation = wksp->run().getGoniometerMatrix();
    rotation.Invert();
  }
  // The "Inelastic" convention is ki-kf and is the default. This matches MDTransfQ3D, and therefore ConvertToMD;
  // note that ConvertToDiffractionMDWorkspace version 1 uses the opposite sign for the same setting.
  const double qSign = (ConfigService::Instance().getString("Q.convention") == "Crystallography") ? -1. : 1.;
  const QConversionTable table(wksp->detectorInfo(), rotation, qSign);

  const double t0 = this->getT0(wksp);
  if (t0 != 0.)
    g_log.information() << "Applying instrument T0 of " << t0 << " microseconds to every event\n";

  H5::H5File h5file(filename, H5F_ACC_RDONLY, Nexus::H5Util::defaultFileAcc());

  std::vector<double> minVals = getProperty(PropertyNames::MIN_VALUES);
  std::vector<double> maxVals = getProperty(PropertyNames::MAX_VALUES);
  if (minVals.empty()) {
    this->progress(.07, "Calculating extents");
    this->calculateExtents(h5file, bankEntryNames, table, t0, minVals, maxVals);
  }
  g_log.information() << "Q extents " << minVals[0] << ":" << maxVals[0] << ", " << minVals[1] << ":" << maxVals[1]
                      << ", " << minVals[2] << ":" << maxVals[2] << "\n";

  // ------------------------------------------------------------------ output workspace
  this->progress(.15, "Creating output workspace");
  auto outputWS = MDEventFactory::CreateMDWorkspace(3, eventType);

  auto frameFactory = Geometry::makeMDFrameFactoryChain();
  const std::string frameName = qSample ? Geometry::QSample::QSampleName : Geometry::QLab::QLabName;
  const std::string dimPrefix = qSample ? "Q_sample_" : "Q_lab_";
  const std::string dimAxes[3] = {"x", "y", "z"};
  for (size_t d = 0; d < 3; ++d) {
    const auto frame = frameFactory->create(Geometry::MDFrameArgument(frameName, ""));
    outputWS->addDimension(std::make_shared<Geometry::MDHistoDimension>(
        dimPrefix + dimAxes[d], "Q" + std::to_string(d + 1), *frame, static_cast<coord_t>(minVals[d]),
        static_cast<coord_t>(maxVals[d]), 1));
  }
  outputWS->setCoordinateSystem(qSample ? Mantid::Kernel::QSample : Mantid::Kernel::QLab);
  outputWS->initialize();

  this->setBoxController(outputWS->getBoxController(), wksp->getInstrument());
  outputWS->splitBox();

  // the instrument, run and sample come along for the ride
  ExperimentInfo_sptr expInfo(wksp->cloneExperimentInfo());
  const std::vector<double> identity{1., 0., 0., 0., 1., 0., 0., 0., 1.};
  expInfo->mutableRun().addProperty("W_MATRIX", identity, true);
  expInfo->mutableRun().addProperty("RUBW_MATRIX", identity, true);
  if (t0 != 0.)
    expInfo->mutableRun().addProperty<double>("T0", t0, true);
  outputWS->addExperimentInfo(expInfo);

  // cache the extents in the workspace's own precision, as the event loop compares against them
  coord_t extentsMin[3];
  coord_t extentsMax[3];
  for (size_t d = 0; d < 3; ++d) {
    extentsMin[d] = outputWS->getDimension(d)->getMinimum();
    extentsMax[d] = outputWS->getDimension(d)->getMaximum();
  }

  // ------------------------------------------------------------------ events
  this->progress(.2, "Reading events");
  if (eventType == EventTypes::FULL) {
    auto ws = std::dynamic_pointer_cast<MDEventWorkspace<MDEvent<3>, 3>>(outputWS);
    this->loadEvents(ws, h5file, bankEntryNames, table, t0, extentsMin, extentsMax);
  } else {
    auto ws = std::dynamic_pointer_cast<MDEventWorkspace<MDLeanEvent<3>, 3>>(outputWS);
    this->loadEvents(ws, h5file, bankEntryNames, table, t0, extentsMin, extentsMax);
  }
  h5file.close();

  if (m_eventsWithNoDetector > 0)
    g_log.warning() << "Discarded " << m_eventsWithNoDetector
                    << " events whose detector id is not in the instrument definition\n";

  setProperty(PropertyNames::OUTPUT_WKSP, outputWS);
}

} // namespace Mantid::DataHandling::LoadEventNexusMD
