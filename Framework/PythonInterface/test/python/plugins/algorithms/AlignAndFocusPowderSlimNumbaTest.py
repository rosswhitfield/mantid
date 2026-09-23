# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2026 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX-License-Identifier: GPL-3.0+
import os
import tempfile
import unittest
from unittest import mock

import h5py
import numpy as np

from mantid.simpleapi import AlignAndFocusPowderSlim, AlignAndFocusPowderSlimNumba, CreateGroupingWorkspace, mtd

VULCAN = "VULCAN_218062.nxs.h5"
L1 = 43.755
L2 = [2.296, 2.070, 2.530]
POLAR = [90.0, 150.0, 65.5]


def write_calibration(filename):
    """A synthetic SaveDiffCal file: three groups, a spread of difc, and some masked pixels.

    It covers a plain range of detector ids, which saves loading the instrument
    to find the real ones; ids that are not in the data change nothing.  The
    range must include every id in the data: AlignAndFocusPowderSlim throws
    when a GroupingWorkspace puts a detector in a group without calibrating it.
    """
    detid = np.arange(0, 600000, dtype=np.int32)
    group = (1 + (detid // 20000) % 3).astype(np.int32)
    # each group's nominal difc, spread by up to 3% so the calibration is not a no-op
    h_over_m = (6.62606896e-34 * 1e10) / (2.0 * 1.674927211e-27 * 1e6)
    nominal = (L1 + np.array(L2)) * np.sin(np.radians(POLAR) / 2.0) / h_over_m
    rng = np.random.default_rng(2024)
    difc = nominal[group - 1] * (1.0 + 0.03 * (rng.random(detid.size) * 2.0 - 1.0))
    use = np.ones(detid.size, dtype=np.int32)
    use[::97] = 0
    with h5py.File(filename, "w") as handle:
        cal = handle.create_group("calibration")
        cal.attrs["NX_class"] = "NXentry"
        for name, values in (("detid", detid), ("difc", difc), ("group", group), ("use", use)):
            cal.create_dataset(name, data=values)
        cal.create_dataset("difa", data=np.zeros(detid.size))
        cal.create_dataset("tzero", data=np.zeros(detid.size))
        # No /calibration/instrument: LoadDiffCal would load VULCAN's and then
        # reject the ids in this range that the instrument does not have.


class AlignAndFocusPowderSlimNumbaTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp_dir = tempfile.TemporaryDirectory()
        cls._cal_file = os.path.join(cls._tmp_dir.name, "vulcan_synthetic_cal.h5")
        write_calibration(cls._cal_file)

        # One group per bank of the current VULCAN definition, then bank 2 unset and
        # bank 3 merged into bank 1. The current definition has banks the file's
        # does not, so some groups name detectors that never appear in the data.
        grouping = CreateGroupingWorkspace(
            InstrumentName="VULCAN", GroupDetectorsBy="bank", OutputWorkspace="vulcan_banks", StoreInADS=False
        )[0]
        for index in range(grouping.getNumberHistograms()):
            group = grouping.readY(index)[0]
            if group == 2:
                grouping.setY(index, np.array([0.0]))
            elif group == 3:
                grouping.setY(index, np.array([1.0]))
        cls._grouping = grouping
        cls._n_groups = len(grouping.getGroupIDs(False))

    @classmethod
    def tearDownClass(cls):
        cls._tmp_dir.cleanup()

    def tearDown(self):
        mtd.clear()

    def _grouped(self):
        """GroupingWorkspace with matching geometry, one L2 and Polar per group."""
        return dict(
            GroupingWorkspace=self._grouping,
            L1=L1,
            L2=list(np.linspace(2.0, 2.5, self._n_groups)),
            Polar=list(np.linspace(60.0, 150.0, self._n_groups)),
        )

    def _compare(self, numba_only=None, calibrated=True, **props):
        """Run both algorithms and require identical histograms and run properties."""
        if "GroupingWorkspace" in props:
            geometry = dict(CalFileName=self._cal_file) if calibrated else {}
        else:
            geometry = dict(CalFileName=self._cal_file, L1=L1, L2=L2, Polar=POLAR) if calibrated else dict(L1=L1, L2=[2.3], Polar=[120.0])
        props = dict(LogAllowList=["frequency", "proton_charge"], **geometry, **props)
        ref = AlignAndFocusPowderSlim(VULCAN, OutputWorkspace="ref", **props)
        new = AlignAndFocusPowderSlimNumba(VULCAN, OutputWorkspace="new", **props, **(numba_only or {}))

        self.assertEqual(new.getNumberHistograms(), ref.getNumberHistograms())
        self.assertGreater(ref.extractY().sum(), 0)
        for index in range(ref.getNumberHistograms()):
            np.testing.assert_array_equal(new.readX(index), ref.readX(index))
            np.testing.assert_array_equal(new.readY(index), ref.readY(index))
            np.testing.assert_array_equal(new.readE(index), ref.readE(index))
        self.assertEqual(new.getAxis(0).getUnit().unitID(), "TOF")

        ref_info, new_info = ref.spectrumInfo(), new.spectrumInfo()
        self.assertAlmostEqual(new_info.l1(), ref_info.l1())
        for index in range(ref.getNumberHistograms()):
            self.assertAlmostEqual(new_info.l2(index), ref_info.l2(index))
            self.assertAlmostEqual(new_info.twoTheta(index), ref_info.twoTheta(index))

        ref_run, new_run = ref.run(), new.run()
        names = sorted(p.name for p in ref_run.getProperties())
        self.assertEqual(sorted(p.name for p in new_run.getProperties()), names)
        for name in names:
            self.assertEqual(str(new_run.getProperty(name).value), str(ref_run.getProperty(name).value), name)
        self.assertEqual(new.getTitle(), ref.getTitle())

    def test_log_dspacing_matches(self):
        self._compare(BinningUnits="dSpacing", BinningMode="Logarithmic", XMin=[0.3], XMax=[3.0], XDelta=[0.0016])

    def test_linear_tof_matches(self):
        self._compare(BinningUnits="TOF", BinningMode="Linear", XMin=[1000.0], XMax=[33000.0], XDelta=[10.0])

    def test_linear_dspacing_matches(self):
        self._compare(BinningUnits="dSpacing", BinningMode="Linear", XMin=[0.3], XMax=[3.0], XDelta=[0.001])

    def test_reads_not_aligned_to_chunks_match(self):
        # an odd read size puts chunk boundaries inside HDF5 chunks, exercising the partial reads
        self._compare(
            numba_only=dict(ReadSizeFromDisk=3000001, ReadAhead=1, ReadThreads=3),
            BinningUnits="dSpacing",
            BinningMode="Logarithmic",
            XMin=[0.3],
            XMax=[3.0],
            XDelta=[0.0016],
        )

    def test_hdf5_fallback_matches(self):
        # files whose event columns are filtered are read through h5py instead
        from plugins.algorithms import AlignAndFocusPowderSlimNumba as module

        with mock.patch.object(module._ChunkMap, "build", return_value=None):
            self._compare(
                numba_only=dict(ReadSizeFromDisk=20000000),
                BinningUnits="dSpacing",
                BinningMode="Logarithmic",
                XMin=[0.3],
                XMax=[3.0],
                XDelta=[0.0016],
            )

    def test_no_calibration_matches(self):
        # difc from the instrument geometry, every detector in one spectrum
        self._compare(calibrated=False, BinningUnits="dSpacing", BinningMode="Logarithmic", XMin=[0.3], XMax=[3.0], XDelta=[0.0016])

    def test_no_calibration_tof_matches(self):
        self._compare(calibrated=False, BinningUnits="TOF", BinningMode="Linear", XMin=[1000.0], XMax=[33000.0], XDelta=[10.0])

    def test_no_calibration_needs_one_spectrum(self):
        with self.assertRaisesRegex(RuntimeError, "one spectrum"):
            AlignAndFocusPowderSlimNumba(VULCAN, OutputWorkspace="new", L1=L1, L2=L2, Polar=POLAR)

    def test_grouping_workspace_overrides_calibration_groups(self):
        self._compare(BinningUnits="dSpacing", BinningMode="Logarithmic", XMin=[0.3], XMax=[3.0], XDelta=[0.0016], **self._grouped())

    def test_grouping_workspace_without_calibration_matches(self):
        self._compare(
            calibrated=False, BinningUnits="dSpacing", BinningMode="Logarithmic", XMin=[0.3], XMax=[3.0], XDelta=[0.0016], **self._grouped()
        )

    def test_grouping_workspace_group_count_must_match_geometry(self):
        with self.assertRaisesRegex(RuntimeError, f"grouping workspace gives {self._n_groups} groups"):
            AlignAndFocusPowderSlimNumba(
                VULCAN, OutputWorkspace="new", CalFileName=self._cal_file, GroupingWorkspace=self._grouping, L1=L1, L2=L2, Polar=POLAR
            )

    def test_ragged_binning_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "ragged"):
            AlignAndFocusPowderSlimNumba(
                VULCAN, OutputWorkspace="new", CalFileName=self._cal_file, L1=L1, L2=L2, Polar=POLAR, XMin=[0.3, 0.4, 0.5]
            )

    def test_group_count_must_match_geometry(self):
        with self.assertRaisesRegex(RuntimeError, "3 groups"):
            AlignAndFocusPowderSlimNumba(VULCAN, OutputWorkspace="new", CalFileName=self._cal_file, L1=L1, L2=[2.0], Polar=[90.0])


if __name__ == "__main__":
    unittest.main()
