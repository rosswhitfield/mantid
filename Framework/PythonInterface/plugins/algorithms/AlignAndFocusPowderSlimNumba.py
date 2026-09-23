# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
"""Prototype of AlignAndFocusPowderSlim in Python, with a numba kernel.

The events never become Mantid objects.  Each bank's ``event_id`` and
``event_time_offset`` columns are read in chunks into two reused NumPy buffers,
histogrammed by a numba kernel into one per-thread array of counts, and only the
final focused histograms become a Workspace2D.  The output gets its geometry from
EditInstrumentGeometry, which is all the C++ algorithm keeps of the instrument,
so the instrument definition is only loaded when there is no calibration file
and difc has to come from the detector positions.

Output is intended to be bit-identical to AlignAndFocusPowderSlim for the
options implemented here.  Floating point order matters for that: the event
TOF is widened to float64 before the multiply, as the C++ loop does, and the
focused difc keeps Mantid's reciprocal form.
"""

import collections
import math
import os
import queue
import threading
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np
from numba import get_num_threads, njit, prange

from mantid.api import (
    AlgorithmFactory,
    FileAction,
    FileProperty,
    MatrixWorkspaceProperty,
    Progress,
    PropertyMode,
    PythonAlgorithm,
    WorkspaceFactory,
)
from mantid.dataobjects import GroupingWorkspaceProperty
from mantid.kernel import (
    Direction,
    FloatArrayBoundedValidator,
    FloatArrayProperty,
    IntBoundedValidator,
    StringArrayProperty,
    StringListValidator,
)

# Kernel::Units::H_OVER_NEUTRON_MASS, from the PhysicalConstants values Mantid compiles in
_H_PLANCK = 6.62606896e-34
_NEUTRON_MASS = 1.674927211e-27
_H_OVER_NEUTRON_MASS = (_H_PLANCK * 1e10) / (2.0 * _NEUTRON_MASS * 1e6)

# createAxisFromRebinParams with full_bins_only=false
_LAST_BIN_COEF = 0.25

_ENTRY = "entry"
_SKIPPED_BANKS = ("bank_error_events", "bank_unmapped_events")
_TIME_UNITS_TO_MICROSECONDS = {
    "microsecond": 1.0,
    "microseconds": 1.0,
    "us": 1.0,
    "second": 1e6,
    "seconds": 1e6,
    "s": 1e6,
    "millisecond": 1e3,
    "milliseconds": 1e3,
    "ms": 1e3,
    "nanosecond": 1e-3,
    "nanoseconds": 1e-3,
    "ns": 1e-3,
}


def _text(value):
    """An HDF5 string attribute or dataset as str, whether stored as bytes, str or an array of them."""
    if isinstance(value, np.ndarray):
        return _text(value.reshape(-1)[0])
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


# ---------------------------------------------------------------------------
# Binning


def _axis_from_rebin_params(xmin, delta, xmax, logarithmic):
    """One range of Kernel::VectorHelper::createAxisFromRebinParams, bit for bit.

    The edges are accumulated, ``x += x * alpha``, not evaluated as
    ``xmin * (1 + alpha) ** k``, and the last edge is ``xmax`` itself.
    """
    alpha = abs(delta)
    edges = [xmin]
    x = xmin
    step = x * alpha if logarithmic else alpha
    while x + (1.0 + _LAST_BIN_COEF) * step <= xmax:
        x = x + step
        edges.append(x)
        step = x * alpha if logarithmic else alpha
    edges.append(xmax)
    return np.array(edges, dtype=np.float64)


def _difc_focused(l1, l2, polar):
    """``1 / tofToDSpacingFactor`` per output spectrum, in Mantid's operation order.

    ``1 / (h / x)`` is not ``x / h`` in floating point; the algebraically equal
    form is an ulp off for some angles, which rescales every bin edge.
    """
    sin_theta = np.sin(np.radians(np.asarray(polar, dtype=np.float64)) / 2.0)
    sin_theta = sin_theta * (l1 + np.asarray(l2, dtype=np.float64))
    return 1.0 / (_H_OVER_NEUTRON_MASS / sin_theta)


def _edges_to_tof(edges, units, difc):
    """What ConvertUnits(Target="TOF") does to one spectrum's bin edges."""
    if units == "TOF":
        return edges.copy()
    if units == "dSpacing":
        return edges * difc
    # MomentumTransfer; the order reverses so the edges stay ascending
    return ((2.0 * math.pi) / edges * difc)[::-1].copy()


# ---------------------------------------------------------------------------
# Kernel


@njit(parallel=True, nogil=True, cache=True)
def _histogram(detid, tof, offset, ratio, spectrum, edges, n_bins, edge0, inv_step, logarithmic, local):
    """Add one chunk of events to ``local``, shaped (threads, spectra, bins).

    Every output spectrum in one pass: ``spectrum`` maps a detector straight to
    its output, so the events are walked once however many spectra there are.
    The bin is computed in closed form and then corrected against the real
    edges, landing where ``upper_bound(edges) - 1`` would.
    """
    n_threads = local.shape[0]
    n = detid.shape[0]
    per_thread = (n + n_threads - 1) // n_threads
    for thread in prange(n_threads):
        start = thread * per_thread
        stop = min(start + per_thread, n)
        for i in range(start, stop):
            index = np.int64(detid[i]) - offset
            if index < 0 or index >= spectrum.shape[0]:
                continue
            spec = spectrum[index]
            if spec < 0:
                continue
            aligned = np.float64(tof[i]) * ratio[index]
            # [front, back), testing the maximum first as ProcessEventsTask does
            if not (aligned < edges[spec, n_bins]):
                continue
            if aligned < edges[spec, 0]:
                continue
            if logarithmic:
                k = np.int64(math.log(aligned / edge0[spec]) * inv_step[spec])
            else:
                k = np.int64((aligned - edge0[spec]) * inv_step[spec])
            if k < 0:
                k = 0
            elif k >= n_bins:
                k = n_bins - 1
            while aligned < edges[spec, k]:
                k -= 1
            while not (aligned < edges[spec, k + 1]):
                k += 1
            local[thread, spec, k] += 1


# ---------------------------------------------------------------------------
# Reading


class _ChunkMap:
    """Where an unfiltered, chunked, 1-D dataset's chunks live in the file."""

    def __init__(self, chunk_elements, itemsize, n_elements, offsets):
        self.chunk_elements = chunk_elements
        self.itemsize = itemsize
        self.n_elements = n_elements
        self.offsets = offsets

    @classmethod
    def build(cls, dataset):
        """The map, or None when the bytes in the file are not the values."""
        dataset_id = dataset.id
        plist = dataset_id.get_create_plist()
        if plist.get_layout() != h5py.h5d.CHUNKED or plist.get_nfilters() != 0:
            return None
        if len(dataset.shape) != 1 or not hasattr(dataset_id, "chunk_iter"):
            return None
        if dataset.dtype.byteorder not in ("=", "|", "<" if np.little_endian else ">"):
            return None
        n_chunks = dataset_id.get_num_chunks()
        if n_chunks == 0:
            return None
        chunk_elements = dataset.chunks[0]
        offsets = np.full(n_chunks, -1, dtype=np.int64)

        def visit(info):
            offsets[info.chunk_offset[0] // chunk_elements] = info.byte_offset

        # chunk_iter walks the index once; get_chunk_info per chunk is quadratic
        dataset_id.chunk_iter(visit)
        if np.any(offsets < 0):
            return None  # unwritten chunks
        return cls(chunk_elements, dataset.dtype.itemsize, dataset.shape[0], offsets)


class _DirectReader:
    """Reads unfiltered event columns straight from the file with ``os.preadv``.

    HDF5 serialises every read behind one library lock and keeps one request
    outstanding; the chunks are stored raw, so reading them ourselves from a
    thread pool gives the disk the queue depth it needs.  Chunks that sit back
    to back in the file are read as one request.
    """

    def __init__(self, filename, threads):
        self._fd = os.open(filename, os.O_RDONLY)
        self._pool = ThreadPoolExecutor(max(1, threads), thread_name_prefix="aafps-io")

    def close(self):
        self._pool.shutdown(wait=True)
        os.close(self._fd)

    def _whole(self, view, offset, nbytes, dest):
        got = os.preadv(self._fd, [view[dest : dest + nbytes]], offset)
        if got != nbytes:
            raise IOError(f"short read: {got} of {nbytes} bytes at {offset}")

    def _part(self, view, offset, nbytes, skip, dest, count):
        raw = os.pread(self._fd, nbytes, offset)
        view[dest : dest + count] = raw[skip : skip + count]

    def submit(self, futures, array, chunk_map, start, stop, dest_element):
        """Queue the reads of elements [start, stop) into ``array[dest_element:]``."""
        view = memoryview(array).cast("B")
        size = chunk_map.itemsize
        chunk_elements = chunk_map.chunk_elements
        chunk_bytes = chunk_elements * size
        offsets = chunk_map.offsets
        run_first = None
        run_last = -1
        run_dest = 0

        def flush():
            if run_first is not None:
                nbytes = (run_last - run_first + 1) * chunk_bytes
                futures.append(self._pool.submit(self._whole, view, int(offsets[run_first]), nbytes, run_dest))

        for chunk in range(start // chunk_elements, (stop - 1) // chunk_elements + 1):
            low_element = chunk * chunk_elements
            low = max(start, low_element)
            high = min(stop, low_element + chunk_elements)
            dest = (dest_element + low - start) * size
            if low == low_element and high == low_element + chunk_elements:
                if run_first is not None and offsets[chunk] == offsets[run_last] + chunk_bytes:
                    run_last = chunk
                else:
                    flush()
                    run_first, run_last, run_dest = chunk, chunk, dest
            else:
                # a partial head or tail chunk, or the padded final one
                flush()
                run_first = None
                futures.append(
                    self._pool.submit(
                        self._part, view, int(offsets[chunk]), chunk_bytes, (low - low_element) * size, dest, (high - low) * size
                    )
                )
        flush()


# error_model="numpy": dividing by zero gives inf, as in C++, rather than raising.
# Monitors on the beam line do that; they are skipped afterwards.
@njit(cache=True, error_model="numpy")
def _difc_from_positions(positions, sample, source, l1):
    """DetectorInfo::difcUncalibrated for every detector, operation for operation.

    ``1 / tofToDSpacingFactor(l1, l2, twoTheta, 0)``, with l2 and twoTheta as
    DetectorInfo computes them from the positions.  Compiled so that sin and
    acos come from the same maths library as the C++, which NumPy's own
    vectorised versions need not match to the last bit.
    """
    beam_x = sample[0] - source[0]
    beam_y = sample[1] - source[1]
    beam_z = sample[2] - source[2]
    beam_norm = math.sqrt(beam_x * beam_x + beam_y * beam_y + beam_z * beam_z)
    difc = np.empty(positions.shape[0])
    for i in range(positions.shape[0]):
        x = positions[i, 0] - sample[0]
        y = positions[i, 1] - sample[1]
        z = positions[i, 2] - sample[2]
        l2 = math.sqrt(x * x + y * y + z * z)
        ratio = (x * beam_x + y * beam_y + z * beam_z) / (l2 * beam_norm)
        if ratio >= 1.0:
            two_theta = 0.0
        elif ratio <= -1.0:
            two_theta = math.pi
        else:
            two_theta = math.acos(ratio)
        sin_theta = math.sin(two_theta / 2)
        sin_theta *= l1 + l2
        difc[i] = 1.0 / ((1.0 * _H_OVER_NEUTRON_MASS) / sin_theta)
    return difc


class _InstrumentDifc:
    """Uncalibrated difc for each detector, and the detectors in each bank.

    Mirrors AlignAndFocusPowderSlim without a calibration: monitors are
    skipped, and a bank's events count only for detectors inside that bank's
    component, which is where getDetectorIDsInBank looks.
    """

    def __init__(self, wksp):
        # componentInfo refers into the workspace, so the workspace must outlive it
        self._wksp = wksp
        detector_info = wksp.detectorInfo()
        self._component_info = wksp.componentInfo()
        self._detids = np.asarray(detector_info.detectorIDs(), dtype=np.int64)
        self._monitor = np.array([detector_info.isMonitor(i) for i in range(len(detector_info))], dtype=bool)
        positions = np.ascontiguousarray(detector_info.allPositions(), dtype=np.float64)
        sample = np.array(list(self._component_info.samplePosition()), dtype=np.float64)
        source = np.array(list(self._component_info.sourcePosition()), dtype=np.float64)
        self._difc = _difc_from_positions(positions, sample, source, detector_info.l1())

    def table(self, bank_name, difc_focus, time_conversion, grouping=None):
        """``(detid_min, ratio, spectrum)`` for one bank, or None when the instrument has no such bank.

        Without ``grouping`` every detector goes to spectrum 0.
        """
        try:
            bank = self._component_info.indexOfAny(bank_name)
        except (RuntimeError, ValueError, IndexError):
            return None
        index = np.asarray(self._component_info.detectorsInSubtree(bank), dtype=np.int64)
        index = index[~self._monitor[index]]
        if index.size == 0:
            return None
        detid = self._detids[index]
        detid_min = int(detid.min())
        ratio = np.zeros(int(detid.max()) - detid_min + 1, dtype=np.float64)
        spectrum = np.full(ratio.size, -1, dtype=np.int32)
        spec = grouping.spectrum_of(detid) if grouping is not None else np.zeros(detid.size, dtype=np.int32)
        keep = spec >= 0
        detid, index, spec = detid[keep], index[keep], spec[keep]
        # difc_focused / difc, as initCalibrationConstants computes it, then the time unit
        ratio[detid - detid_min] = (difc_focus[spec] / self._difc[index]) * time_conversion
        spectrum[detid - detid_min] = spec
        return detid_min, ratio, spectrum


class _Grouping:
    """A GroupingWorkspace as a detector id to output spectrum lookup.

    Spectra are the non-zero group ids in ascending order, as
    ``getGroupIDs(false)`` gives them to AlignAndFocusPowderSlim.  A detector
    listed in more than one group lands in the last, as the C++ map does when
    it is filled group by group.
    """

    def __init__(self, wksp):
        self.group_ids = np.asarray(wksp.getGroupIDs(False), dtype=np.int64)
        if self.group_ids.size == 0:
            raise RuntimeError("The grouping workspace has no groups")
        members = [np.asarray(wksp.getDetectorIDsOfGroup(int(group)), dtype=np.int64) for group in self.group_ids]
        everything = np.concatenate(members)
        if everything.size == 0:
            raise RuntimeError("The grouping workspace has no detectors in its groups")
        self._offset = int(everything.min())
        self._lookup = np.full(int(everything.max()) - self._offset + 1, -1, dtype=np.int32)
        for spectrum, detids in enumerate(members):
            self._lookup[detids - self._offset] = spectrum

    @property
    def n_spectra(self):
        return int(self.group_ids.size)

    def spectrum_of(self, detids):
        """The output spectrum of each detector id, or -1 for one in no group."""
        index = np.asarray(detids, dtype=np.int64) - self._offset
        inside = (index >= 0) & (index < self._lookup.size)
        spectrum = np.full(index.size, -1, dtype=np.int32)
        spectrum[inside] = self._lookup[index[inside]]
        return spectrum


class _Buffers:
    """One chunk's detector ids and times of flight, reused for every chunk."""

    def __init__(self, capacity):
        self.detid = np.empty(capacity, dtype=np.uint32)
        self.tof = np.empty(capacity, dtype=np.float32)
        self.size = 0


# ---------------------------------------------------------------------------
# The algorithm


class AlignAndFocusPowderSlimNumba(PythonAlgorithm):
    def category(self):
        return "Diffraction\\Reduction"

    def seeAlso(self):
        return ["AlignAndFocusPowderSlim"]

    def summary(self):
        return "Prototype: AlignAndFocusPowderSlim in Python with a numba kernel, for comparing speed."

    def PyInit(self):
        self.declareProperty(FileProperty("Filename", "", action=FileAction.Load, extensions=[".nxs.h5"]), doc="Raw event NeXus file")
        self.declareProperty(
            FileProperty("CalFileName", "", action=FileAction.OptionalLoad, extensions=[".h5"]),
            doc="Calibration file from SaveDiffCal, giving difc, grouping and mask. Without one, difc comes from the "
            "instrument geometry and every detector is focused into one spectrum.",
        )
        self.declareProperty(
            GroupingWorkspaceProperty("GroupingWorkspace", "", direction=Direction.Input, optional=PropertyMode.Optional),
            doc="Grouping of detectors into output spectra. Takes precedence over the grouping in CalFileName.",
        )
        positive = FloatArrayBoundedValidator(lower=0.0, exclusive=True)
        self.declareProperty(FloatArrayProperty("XMin", [0.1], positive), doc="Minimum x-value of the output binning")
        self.declareProperty(FloatArrayProperty("XDelta", [0.0016]), doc="Bin size of the output binning")
        self.declareProperty(FloatArrayProperty("XMax", [2.0], positive), doc="Maximum x-value of the output binning")
        self.declareProperty(
            "BinningUnits", "dSpacing", StringListValidator(["dSpacing", "TOF", "MomentumTransfer"]), doc="Units of XMin, XDelta and XMax"
        )
        self.declareProperty("BinningMode", "Logarithmic", StringListValidator(["Logarithmic", "Linear"]), doc="Binning behaviour")
        self.declareProperty("L1", 0.0, doc="Primary flight path of the output instrument")
        self.declareProperty(FloatArrayProperty("L2", []), doc="Secondary flight path of each output spectrum")
        self.declareProperty(FloatArrayProperty("Polar", []), doc="Polar angle (two-theta) of each output spectrum, in degrees")
        self.declareProperty(FloatArrayProperty("Azimuthal", []), doc="Azimuthal angle of each output spectrum, in degrees")
        self.declareProperty(StringArrayProperty("LogAllowList", []), doc="If given, load only these logs")
        self.declareProperty(
            StringArrayProperty("LogBlockList", ["Phase\\*", "Speed\\*", "BL\\*:Chop:\\*", "chopper\\*TDC"]),
            doc="Logs not to load; ignored when LogAllowList is given",
        )
        self.declareProperty("ReadSizeFromDisk", 8000000, IntBoundedValidator(lower=1), doc="Number of events read from disk at a time")
        self.declareProperty(
            "ReadAhead", 4, IntBoundedValidator(lower=1), doc="Chunks whose reads may be in flight while one is histogrammed"
        )
        self.declareProperty("ReadThreads", 16, IntBoundedValidator(lower=1), doc="Threads issuing reads for the direct reader")
        self.declareProperty(
            MatrixWorkspaceProperty("OutputWorkspace", "", direction=Direction.Output), doc="Focused workspace, in time-of-flight"
        )

    def validateInputs(self):
        issues = {}
        for name in ("XMin", "XDelta", "XMax"):
            if len(self.getProperty(name).value) != 1:
                issues[name] = "The prototype supports one value only, not ragged binning"
        delta = self.getProperty("XDelta").value
        if len(delta) == 1 and (delta[0] == 0.0 or not math.isfinite(delta[0])):
            issues["XDelta"] = "XDelta must be finite and non-zero"
        if len(self.getProperty("XMin").value) == 1 and len(self.getProperty("XMax").value) == 1:
            if not self.getProperty("XMax").value[0] > self.getProperty("XMin").value[0]:
                issues["XMax"] = "XMax must exceed XMin"
        if not self.getProperty("L1").value > 0.0:
            issues["L1"] = "L1 must be positive"
        n_l2 = len(self.getProperty("L2").value)
        if n_l2 == 0 or n_l2 != len(self.getProperty("Polar").value):
            issues["Polar"] = "L2 and Polar must be given, one value per output spectrum"
        azimuthal = self.getProperty("Azimuthal").value
        if len(azimuthal) not in (0, n_l2):
            issues["Azimuthal"] = "Azimuthal must be empty or match L2"
        if len(self.getProperty("LogAllowList").value) and not self.getProperty("LogBlockList").isDefault:
            issues["LogAllowList"] = "Cannot specify both allow and block lists"
        return issues

    def PyExec(self):
        filename = self.getProperty("Filename").value
        l1 = self.getProperty("L1").value
        l2 = self.getProperty("L2").value
        polar = self.getProperty("Polar").value
        units = self.getProperty("BinningUnits").value
        logarithmic = self.getProperty("BinningMode").value == "Logarithmic"
        read_size = self.getProperty("ReadSizeFromDisk").value

        difc_focus = _difc_focused(l1, l2, polar)
        grouping_wksp = self.getProperty("GroupingWorkspace").value
        grouping = _Grouping(grouping_wksp) if grouping_wksp is not None else None
        cal_filename = self.getProperty("CalFileName").value
        if cal_filename:
            detid_min, spectrum_of_row, difc, detid, n_spectra = self._load_calibration(cal_filename, grouping)
            source = "grouping workspace" if grouping else "calibration file"
            if n_spectra != len(l2):
                raise RuntimeError(f"The {source} gives {n_spectra} groups but L2 and Polar have {len(l2)} values")

            def detector_table(bank_name, time_conversion):
                return (detid_min, *self._detector_table(detid, detid_min, spectrum_of_row, difc, difc_focus, time_conversion))

        else:
            n_spectra = grouping.n_spectra if grouping else 1
            if n_spectra != len(l2):
                if grouping:
                    raise RuntimeError(f"The grouping workspace gives {n_spectra} groups but L2 and Polar have {len(l2)} values")
                raise RuntimeError(
                    "Without a calibration file or grouping workspace every detector is focused into one spectrum, "
                    f"but L2 and Polar have {len(l2)} values"
                )
            instrument = _InstrumentDifc(self._load_instrument(filename))

            def detector_table(bank_name, time_conversion):
                return instrument.table(bank_name, difc_focus, time_conversion, grouping)

        delta = self.getProperty("XDelta").value[0]
        d_edges = _axis_from_rebin_params(
            self.getProperty("XMin").value[0], delta if not logarithmic else -abs(delta), self.getProperty("XMax").value[0], logarithmic
        )
        edges = np.stack([_edges_to_tof(d_edges, units, difc_focus[s]) for s in range(n_spectra)])
        n_bins = edges.shape[1] - 1
        # The closed-form guess of the bin assumes the edges are uniform in the
        # value or its logarithm, which scaling by difc preserves.  It is only a
        # guess: the kernel corrects it against the edges, so a poor one costs
        # time, not correctness.
        if units == "MomentumTransfer":
            inv_step = np.zeros(n_spectra)  # reversed edges have no simple form
        elif logarithmic:
            inv_step = np.full(n_spectra, 1.0 / math.log1p(abs(delta)))
        elif units == "dSpacing":
            inv_step = 1.0 / (abs(delta) * difc_focus)
        else:
            inv_step = np.full(n_spectra, 1.0 / abs(delta))

        with h5py.File(filename, "r") as handle:
            banks = self._banks(handle)
            n_threads = get_num_threads()
            local = np.zeros((n_threads, n_spectra, n_bins), dtype=np.int64)
            reader = _DirectReader(filename, self.getProperty("ReadThreads").value)
            try:
                progress = Progress(self, start=0.0, end=0.95, nreports=len(banks))
                for entry_name, n_events, time_conversion in banks:
                    progress.report(f"Processing {entry_name}")
                    table = detector_table(entry_name.split("_", 1)[0], time_conversion)
                    if table is None:
                        self.log().warning(f"Skipping {entry_name}: the instrument has no component with that bank's name")
                        continue
                    bank_detid_min, ratio, spectrum = table
                    group = handle[f"{_ENTRY}/{entry_name}"]
                    self._stream(
                        group, n_events, read_size, reader, bank_detid_min, ratio, spectrum, edges, n_bins, inv_step, logarithmic, local
                    )
            finally:
                reader.close()

        counts = local.sum(axis=0)
        self.setProperty("OutputWorkspace", self._make_workspace(filename, edges, counts, l1, l2, polar))

    # -----------------------------------------------------------------------

    def _load_calibration(self, filename, grouping=None):
        """``(detid_min, spectrum_per_row, difc_per_row, detid_per_row, n_spectra)`` from a SaveDiffCal file.

        The output spectrum comes from ``grouping`` when given, else it is the
        rank of the row's group id among the non-zero ids, as the grouping map
        in AlignAndFocusPowderSlim.  Rows in no group, with ``use == 0``, or
        with zero difc are ignored (spectrum -1).
        """
        with h5py.File(filename, "r") as handle:
            cal = handle["calibration"]
            detid = np.asarray(cal["detid"][:], dtype=np.int64)
            difc = np.asarray(cal["difc"][:], dtype=np.float64)
            group = np.asarray(cal["group"][:], dtype=np.int64) if "group" in cal else np.ones(detid.size, dtype=np.int64)
            use = np.asarray(cal["use"][:], dtype=np.int64) if "use" in cal else np.ones(detid.size, dtype=np.int64)
        if grouping is not None:
            spectrum = grouping.spectrum_of(detid)
            n_spectra = grouping.n_spectra
        else:
            group_ids = np.unique(group[group != 0])
            spectrum = np.searchsorted(group_ids, group).astype(np.int32)
            spectrum[group == 0] = -1
            n_spectra = group_ids.size
        spectrum[(use == 0) | (difc == 0.0)] = -1
        return int(detid.min()), spectrum, difc, detid, n_spectra

    def _load_instrument(self, filename):
        """A one-spectrum workspace carrying the instrument, loaded as LoadEventNexus::loadInstrument does.

        Only needed without a calibration file, for the detector positions.
        """
        wksp = WorkspaceFactory.create("Workspace2D", NVectors=1, XLength=2, YLength=1)
        from_nexus = self.createChildAlgorithm("LoadIDFFromNexus", enableLogging=False)
        from_nexus.setProperty("Workspace", wksp)
        from_nexus.setProperty("Filename", filename)
        from_nexus.setProperty("InstrumentParentPath", _ENTRY)
        try:
            from_nexus.execute()
            wksp = from_nexus.getProperty("Workspace").value
        except RuntimeError:
            # no instrument_xml in the file: find the definition by name and run start
            with h5py.File(filename, "r") as handle:
                name = _text(handle[f"{_ENTRY}/instrument/name"][()])
                if "start_time" in handle[_ENTRY]:
                    wksp.mutableRun().addProperty("run_start", _text(handle[f"{_ENTRY}/start_time"][()]), True)
            by_name = self.createChildAlgorithm("LoadInstrument")
            by_name.setProperty("Workspace", wksp)
            by_name.setProperty("InstrumentName", name)
            by_name.setProperty("RewriteSpectraMap", False)
            by_name.execute()
            wksp = by_name.getProperty("Workspace").value
        return wksp

    def _detector_table(self, detid, detid_min, spectrum_of_row, difc, difc_focus, time_conversion):
        """Dense lookups indexed by ``detid - detid_min``: the calibration factor and output spectrum."""
        size = int(detid.max()) - detid_min + 1
        ratio = np.zeros(size, dtype=np.float64)
        spectrum = np.full(size, -1, dtype=np.int32)
        keep = spectrum_of_row >= 0
        index = detid[keep] - detid_min
        spec = spectrum_of_row[keep]
        spectrum[index] = spec
        # difc_focused / difc, as initCalibrationConstantsFromCalWS computes it,
        # then scaled by the time unit as BankCalibration does
        ratio[index] = (difc_focus[spec] / difc[keep]) * time_conversion
        return ratio, spectrum

    def _banks(self, handle):
        """``(entry_name, n_events, time_conversion)`` for each NXevent_data under the entry."""
        banks = []
        for name, group in handle[_ENTRY].items():
            if not isinstance(group, h5py.Group) or name in _SKIPPED_BANKS:
                continue
            if _text(group.attrs.get("NX_class", b"")) != "NXevent_data":
                continue
            if "event_time_offset" not in group:
                continue
            tof = group["event_time_offset"]
            units = _text(tof.attrs.get("units", b"microsecond")).lower()
            if units not in _TIME_UNITS_TO_MICROSECONDS:
                raise RuntimeError(f"{name}/event_time_offset has unsupported units {units!r}")
            banks.append((name, int(tof.shape[0]), _TIME_UNITS_TO_MICROSECONDS[units]))
        if not banks:
            raise RuntimeError("No NXevent_data entries found in file")
        return banks

    def _stream(self, group, n_events, read_size, reader, detid_min, ratio, spectrum, edges, n_bins, inv_step, logarithmic, local):
        """Read the bank a chunk at a time and histogram it, reading ahead.

        A reader thread keeps up to ``ReadAhead`` chunks of reads in flight
        while the kernel, which releases the GIL, bins the chunk before them.
        The buffers are reused, so memory is bounded by ``ReadAhead + 1`` chunks.
        """
        if n_events == 0:
            return
        detid_dataset = group["event_id"]
        tof_dataset = group["event_time_offset"]
        detid_map = _ChunkMap.build(detid_dataset)
        tof_map = _ChunkMap.build(tof_dataset)
        # raw bytes go straight into the buffers, so the file types must match them
        direct = (
            detid_map is not None
            and tof_map is not None
            and detid_dataset.dtype == np.dtype(np.uint32)
            and tof_dataset.dtype == np.dtype(np.float32)
        )
        capacity = min(read_size, n_events)
        read_ahead = self.getProperty("ReadAhead").value
        free = queue.Queue()
        pending = queue.Queue()
        # one buffer being histogrammed, the rest being read
        for _ in range(read_ahead + 1):
            free.put(_Buffers(capacity))
        edge0 = np.ascontiguousarray(edges[:, 0])

        stop_reading = threading.Event()

        def fill():
            in_flight = collections.deque()

            def complete():
                buffers, futures = in_flight.popleft()
                for future in futures:
                    future.result()
                pending.put(buffers)

            try:
                for start in range(0, n_events, read_size):
                    stop = min(start + read_size, n_events)
                    # hand over finished chunks rather than wait for a buffer they hold
                    while in_flight and free.empty():
                        complete()
                    buffers = None
                    while buffers is None:
                        if stop_reading.is_set():
                            return  # the consumer failed and no longer wants chunks
                        try:
                            buffers = free.get(timeout=0.1)
                        except queue.Empty:
                            pass
                    buffers.size = stop - start
                    futures = []
                    if direct:
                        reader.submit(futures, buffers.detid, detid_map, start, stop, 0)
                        reader.submit(futures, buffers.tof, tof_map, start, stop, 0)
                    else:
                        # filtered, contiguous or other types: let HDF5 decode and convert
                        detid_dataset.read_direct(buffers.detid, np.s_[start:stop], np.s_[0 : buffers.size])
                        tof_dataset.read_direct(buffers.tof, np.s_[start:stop], np.s_[0 : buffers.size])
                    in_flight.append((buffers, futures))
                while in_flight:
                    complete()
            except BaseException as error:
                pending.put(error)
            else:
                pending.put(None)

        thread = threading.Thread(target=fill, name="aafps-reader", daemon=True)
        thread.start()
        try:
            while True:
                item = pending.get()
                if item is None:
                    break
                if isinstance(item, BaseException):
                    raise item
                try:
                    _histogram(
                        item.detid[: item.size],
                        item.tof[: item.size],
                        np.int64(detid_min),
                        ratio,
                        spectrum,
                        edges,
                        n_bins,
                        edge0,
                        inv_step,
                        logarithmic,
                        local,
                    )
                finally:
                    free.put(item)
        finally:
            # if the kernel raised, the reader may be waiting for a buffer
            stop_reading.set()
            thread.join()

    def _load_entry_metadata(self, filename, wksp):
        """What LoadEventNexus::loadEntryMetadata puts on the workspace."""
        with h5py.File(filename, "r") as handle:
            entry = handle[_ENTRY]

            def string(name):
                if name in entry and entry[name].dtype.kind in "SOU":
                    return _text(entry[name][()])
                return ""

            if title := string("title"):
                wksp.setTitle(title)
            run = wksp.mutableRun()
            if notes := string("notes"):
                run.addProperty("file_notes", notes, True)
            run_number = string("run_number")
            if not run_number and "run_number" in entry and entry["run_number"].dtype.kind in "iu":
                run_number = str(int(np.asarray(entry["run_number"][()]).reshape(-1)[0]))
            if run_number:
                run.addProperty("run_number", run_number, True)
            if experiment := string("experiment_identifier"):
                run.addProperty("experiment_identifier", experiment, True)
            # The sample name is not set: Python has no setter for it.
            if "duration" in entry:
                duration = np.asarray(entry["duration"][()], dtype=np.float64).reshape(-1)
                if duration.size == 1:
                    units = _text(entry["duration"].attrs.get("units", b""))
                    run.addProperty("duration", float(duration[0]), units, True)

    def _make_workspace(self, filename, edges, counts, l1, l2, polar):
        """The only Mantid workspace this algorithm makes: the focused histograms."""
        n_spectra, n_edges = edges.shape
        wksp = WorkspaceFactory.create("Workspace2D", NVectors=n_spectra, XLength=n_edges, YLength=n_edges - 1)
        for index in range(n_spectra):
            wksp.setX(index, edges[index])
            wksp.setY(index, counts[index].astype(np.float64))
            wksp.setE(index, np.sqrt(counts[index].astype(np.float64)))
            wksp.getSpectrum(index).setSpectrumNo(index + 1)
        wksp.getAxis(0).setUnit("TOF")
        wksp.setYUnit("Counts")
        wksp.setDistribution(False)

        self._load_entry_metadata(filename, wksp)
        if self.getProperty("BinningUnits").value != "TOF":
            # ConvertUnits records the mode when AlignAndFocusPowderSlim converts its bin edges
            wksp.mutableRun().addProperty("deltaE-mode", "Elastic", True)

        logs = self.createChildAlgorithm("LoadNexusLogs", enableLogging=False)
        logs.setProperty("Workspace", wksp)
        logs.setProperty("Filename", filename)
        allow = self.getProperty("LogAllowList").value
        if len(allow):
            logs.setProperty("AllowList", list(allow))
        else:
            logs.setProperty("BlockList", list(self.getProperty("LogBlockList").value))
        logs.execute()
        # LoadEventNexus::runLoadNexusLogs then records when the run started
        run = wksp.mutableRun()
        if run.hasProperty("proton_charge") and run.getProperty("proton_charge").size() > 0:
            run.addProperty("run_start", run.getProperty("proton_charge").firstTime().toISO8601String(), True)
        elif run.hasProperty("start_time"):
            run.addProperty("run_start", run.getProperty("start_time").value, True)

        with h5py.File(filename, "r") as handle:
            name = _text(handle[f"{_ENTRY}/instrument/name"][()])
        edit = self.createChildAlgorithm("EditInstrumentGeometry", enableLogging=False)
        edit.setProperty("Workspace", wksp)
        edit.setProperty("PrimaryFlightPath", l1)
        edit.setProperty("L2", list(l2))
        edit.setProperty("Polar", list(polar))
        azimuthal = self.getProperty("Azimuthal").value
        edit.setProperty("Azimuthal", list(azimuthal) if len(azimuthal) else [0.0] * len(l2))
        edit.setProperty("InstrumentName", name)
        edit.execute()
        return edit.getProperty("Workspace").value


AlgorithmFactory.subscribe(AlignAndFocusPowderSlimNumba)
