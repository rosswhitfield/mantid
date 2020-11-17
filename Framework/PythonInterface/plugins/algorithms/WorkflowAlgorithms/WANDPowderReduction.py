# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
from mantid.api import (
    ADSValidator,
    AlgorithmFactory,
    DataProcessorAlgorithm,
    MatrixWorkspaceProperty,
    PropertyMode,
    IEventWorkspace,
)
from mantid.dataobjects import MaskWorkspaceProperty
from mantid.simpleapi import (
    AnalysisDataService,
    ConvertSpectrumAxis,
    ConjoinWorkspaces,
    Transpose,
    ResampleX,
    CopyInstrumentParameters,
    CloneWorkspace,
    Divide,
    DeleteWorkspace,
    Scale,
    MaskAngle,
    ExtractMask,
    Minus,
    SumSpectra,
    ExtractUnmaskedSpectra,
    mtd,
    BinaryOperateMasks,
    Integration,
)
from mantid.kernel import (
    StringListValidator,
    Direction,
    Property,
    FloatBoundedValidator,
    StringArrayProperty,
)


class WANDPowderReduction(DataProcessorAlgorithm):
    temp_workspace_list = [
        "__ws_conjoined",
        "_ws_cal",
        "_ws_cal_resampled",
        "_ws_tmp",
        "_ws_resampled",
    ]

    def category(self):
        return "Diffraction\\Reduction"

    def seeAlso(self):
        return ["LoadWAND", "SaveFocusedXYE"]

    def summary(self):
        return "Performs powder diffraction data reduction for WAND"

    def name(self):
        return "WANDPowderReduction"

    def validateInputs(self):
        issues = dict()
        return issues

    def PyInit(self):
        self.declareProperty(
            StringArrayProperty(
                "InputWorkspace", direction=Direction.Input, validator=ADSValidator(),
            ),
            doc="The main input workspace[s].",
        )

        self.declareProperty(
            "BackgroundWorkspace",
            "",
            direction=Direction.Input,
            doc="The background workspace to be subtracted.",
        )

        self.declareProperty(
            MatrixWorkspaceProperty(
                "CalibrationWorkspace",
                "",
                optional=PropertyMode.Optional,
                direction=Direction.Input,
            ),
            doc="The calibration (vandiaum) workspace.",
        )

        self.declareProperty(
            "BackgroundScale",
            1.0,
            validator=FloatBoundedValidator(0.0),
            doc="The background will be scaled by this number before being subtracted.",
        )

        self.declareProperty(
            MaskWorkspaceProperty(
                "MaskWorkspace",
                "",
                optional=PropertyMode.Optional,
                direction=Direction.Input,
            ),
            doc="The mask from this workspace will be applied before reduction",
        )

        self.copyProperties("ConvertSpectrumAxis", ["Target", "EFixed"])

        self.copyProperties("ResampleX", ["XMin", "XMax", "NumberBins", "LogBinning"])

        self.declareProperty(
            "NormaliseBy",
            "Monitor",
            StringListValidator(["None", "Time", "Monitor"]),
            "Normalise to monitor or time. ",
        )

        self.declareProperty(
            "MaskAngle",
            Property.EMPTY_DBL,
            "Phi angle above which will be masked. See :ref:`MaskAngle <algm-MaskAngle>` for details.",
        )

        self.declareProperty(
            MatrixWorkspaceProperty("OutputWorkspace", "", direction=Direction.Output),
            doc="Output Workspace",
        )

    def PyExec(self):
        data = self.getProperty("InputWorkspace").value  # [1~n]
        bkg = self.getProperty("BackgroundWorkspace").value  # [1~n]
        cal = self.getProperty("CalibrationWorkspace").value  # [1]
        numberBins = self.getProperty("NumberBins").value
        outWS = self.getPropertyValue("OutputWorkspace")

        # NOTE:
        # StringArrayProperty cannot be optional, so the background can only be passed in as a string
        # or a list, which will be manually unpacked here
        if bkg != "":
            bkg = [
                AnalysisDataService.retrieve(me)
                for me in map(str.strip, bkg.split(","))
            ]

        # convert all of the input workspaces into spectrum of "target" units (generally angle)
        data, masks = self._convert_data(data)

        # determine x-range
        xMin, xMax = self._locate_global_xlimit(data)

        # BEGIN_FOR: prcess_spectra
        for n, (_wsn, _mskn) in enumerate(zip(data, masks)):
            # resample spectra
            _ws_resampled = ResampleX(
                InputWorkspace=_wsn,
                XMin=xMin,
                XMax=xMax,
                NumberBins=numberBins,
                EnableLogging=False,
            )

            # calibration
            if cal is not None:
                _ws_cal_resampled = self._resample_calibration(_wsn, _mskn, xMin, xMax)
                _ws_resampled = Divide(
                    LHSWorkspace=_ws_resampled,
                    RHSWorkspace=_ws_cal_resampled,
                    EnableLogging=False,
                )
            else:
                _ws_cal_resampled = None

            _ws_resampled = Scale(
                InputWorkspace=_ws_resampled,
                Factor=self._get_scale(cal) / self._get_scale(_wsn),
                EnableLogging=False,
            )

            # background
            if bkg != "":
                bgn = bkg[n] if isinstance(bkg, list) else bkg

                _ws_bkg_resampled = self._resample_background(
                    bgn, _wsn, _mskn, xMin, xMax, _ws_cal_resampled
                )

                _ws_resampled = Minus(
                    LHSWorkspace=_ws_resampled,
                    RHSWorkspace=_ws_bkg_resampled,
                    EnableLogging=False,
                )

            # conjoin
            if n < 1:
                CloneWorkspace(
                    InputWorkspace=_ws_resampled,
                    OutputWorkspace="__ws_conjoined",
                    EnableLogging=False,
                )
            else:
                ConjoinWorkspaces(
                    InputWorkspace1="__ws_conjoined",
                    InputWorkspace2=_ws_resampled,
                    CheckOverlapping=False,
                    EnableLogging=False,
                )
        # END_FOR: prcess_spectra

        # Step_3: sum all spectra
        # ref: https://docs.mantidproject.org/nightly/algorithms/SumSpectra-v1.html
        if cal is not None:
            SumSpectra(
                InputWorkspace="__ws_conjoined",
                OutputWorkspace=outWS,
                WeightedSum=True,
                MultiplyBySpectra=False,
                EnableLogging=False,
            )
        else:
            SumSpectra(
                InputWorkspace="__ws_conjoined",
                OutputWorkspace=outWS,
                WeightedSum=True,
                MultiplyBySpectra=True,
                EnableLogging=False,
            )

        self.setProperty("OutputWorkspace", outWS)

        # Step_4: remove temp workspaces
        [
            DeleteWorkspace(ws, EnableLogging=False)
            for ws in self.temp_workspace_list
            if mtd.doesExist(ws)
        ]

    def _get_scale(self, x):
        """return the scale factor needed during normalization"""
        normaliseBy = self.getProperty("NormaliseBy").value

        if x is None:
            return 1
        else:
            if str(normaliseBy).lower() == "none":
                return 1
            elif str(normaliseBy).lower() == "monitor":
                return mtd[str(x)].run().getProtonCharge()
            elif str(normaliseBy).lower() == "time":
                return mtd[str(x)].run().getLogData("duration").value
            else:
                raise ValueError(f"Unknown normalize type: {normaliseBy}")

    def _convert_data(self, input_workspaces):
        mask = self.getProperty("MaskWorkspace").value
        mask_angle = self.getProperty("MaskAngle").value
        target = self.getProperty("Target").value
        e_fixed = self.getProperty("EFixed").value

        # NOTE:
        # Due to range difference among incoming spectra, a common bin para is needed
        # such that all data can be binned exactly the same way.

        # BEGIN_FOR: located_global_xMin&xMax
        output_workspaces = []  # names will
        mask_workspaces = []
        for n, _wksp_in in enumerate(input_workspaces):
            _mask_n = f'__mask_{n}'  # mask for n-th
            _wksp_out = f'__wskp_{n}'    # output workspace for n-th
            self.temp_workspace_list.append(_mask_n)  # cleanup later

            ExtractMask(InputWorkspace=_wksp_in, OutputWorkspace=_mask_n, EnableLogging=False)
            if mask_angle != Property.EMPTY_DBL:
                MaskAngle(
                    Workspace=_mask_n,
                    MinAngle=mask_angle,
                    Angle="Phi",
                    EnableLogging=False,
                )
            if mask is not None:
                # might be a bug if the mask angle isn't set
                BinaryOperateMasks(
                    InputWorkspace1=_mask_n,
                    InputWorkspace2=mask,
                    OperationType='OR',
                    OutputWorkspace=_mask_n,
                    EnableLogging=False,
                )

            ExtractUnmaskedSpectra(
                InputWorkspace=_wksp_in, OutputWorkspace=_wksp_out,
                MaskWorkspace=_mask_n, EnableLogging=False
            )
            if isinstance(mtd[_wksp_out], IEventWorkspace):
                Integration(InputWorkspace=_wksp_out, OutputWorkspace=_wksp_out, EnableLogging=False)
            ConvertSpectrumAxis(
                InputWorkspace=_wksp_out,
                OutputWorkspace=_wksp_out,
                Target=target,
                EFixed=e_fixed,
                EnableLogging=False,
            )
            Transpose(
                InputWorkspace=_wksp_out, OutputWorkspace=_wksp_out, EnableLogging=False
            )

            # append to the list of processed workspaces
            output_workspaces.append(_wksp_out)
            mask_workspaces.append(_mask_n)

        return output_workspaces, mask_workspaces

    def _locate_global_xlimit(self, workspaces):
        """Find the global bin from all spectrum"""
        # Due to range difference among incoming spectra, a common bin para is needed
        # such that all data can be binned exactly the same way.

        # use the supplied start value
        if self.getProperty('xMin').isDefault:
            _xMin = 1e16
        else:
            _xMin = self.getProperty("XMin").value
        if self.getProperty('xMax').isDefault:
            _xMax = -1e16
        else:
            _xMax = self.getProperty("XMax").value
        # if both were set there is nothing to do
        if _xMin < _xMax and _xMin < 1e16 and _xMax > -1e16:
            return _xMin, _xMax

        # update values based on all workspaces
        for name in workspaces:
            _ws_tmp = mtd[name]
            _xMin = min(_xMin, _ws_tmp.readX(0).min())
            _xMax = max(_xMax, _ws_tmp.readX(0).max())

        return _xMin, _xMax

    def _resample_calibration(
        self, current_workspace, mask_name, x_min, x_max,
    ):
        """Perform resample on calibration"""
        cal = self.getProperty("CalibrationWorkspace").value
        target = self.getProperty("Target").value
        e_fixed = self.getProperty("EFixed").value
        number_bins = self.getProperty("NumberBins").value

        _ws_cal = ExtractUnmaskedSpectra(
            InputWorkspace=cal, MaskWorkspace=mask_name, EnableLogging=False
        )
        if isinstance(mtd["_ws_cal"], IEventWorkspace):
            _ws_cal = Integration(InputWorkspace=_ws_cal, EnableLogging=False)
        CopyInstrumentParameters(
            InputWorkspace=current_workspace,
            OutputWorkspace=_ws_cal,
            EnableLogging=False,
        )
        _ws_cal = ConvertSpectrumAxis(
            InputWorkspace=_ws_cal, Target=target, EFixed=e_fixed, EnableLogging=False,
        )
        _ws_cal = Transpose(InputWorkspace=_ws_cal, EnableLogging=False)
        return ResampleX(
            InputWorkspace=_ws_cal,
            XMin=x_min,
            XMax=x_max,
            NumberBins=number_bins,
            EnableLogging=False,
        )

    def _resample_background(
        self,
        current_background,
        current_workspace,
        make_name,
        x_min,
        x_max,
        resmapled_calibration,
    ):
        """Perform resample on given background"""
        cal = self.getProperty("CalibrationWorkspace").value
        target = self.getProperty("Target").value
        e_fixed = self.getProperty("EFixed").value
        number_bins = self.getProperty("NumberBins").value

        _ws_bkg = ExtractUnmaskedSpectra(
            InputWorkspace=current_background,
            MaskWorkspace=make_name,
            EnableLogging=False,
        )

        if isinstance(mtd["_ws_bkg"], IEventWorkspace):
            _ws_bkg = Integration(InputWorkspace=_ws_bkg, EnableLogging=False)

        CopyInstrumentParameters(
            InputWorkspace=current_workspace,
            OutputWorkspace=_ws_bkg,
            EnableLogging=False,
        )

        _ws_bkg = ConvertSpectrumAxis(
            InputWorkspace=_ws_bkg, Target=target, EFixed=e_fixed, EnableLogging=False,
        )

        _ws_bkg = Transpose(InputWorkspace=_ws_bkg, EnableLogging=False)

        _ws_bkg_resampled = ResampleX(
            InputWorkspace=_ws_bkg,
            XMin=x_min,
            XMax=x_max,
            NumberBins=number_bins,
            EnableLogging=False,
        )

        if cal is not None:
            _ws_bkg_resampled = Divide(
                LHSWorkspace=_ws_bkg_resampled,
                RHSWorkspace=resmapled_calibration,
                EnableLogging=False,
            )

        _ws_bkg_resampled = Scale(
            InputWorkspace=_ws_bkg_resampled,
            Factor=self._get_scale(cal) / self._get_scale(current_background),
            EnableLogging=False,
        )

        _ws_bkg_resampled = Scale(
            InputWorkspace=_ws_bkg_resampled,
            Factor=self.getProperty("BackgroundScale").value,
            EnableLogging=False,
        )

        return _ws_bkg_resampled


AlgorithmFactory.subscribe(WANDPowderReduction)
