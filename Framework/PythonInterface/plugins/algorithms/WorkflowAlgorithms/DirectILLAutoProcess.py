# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

import DirectILL_common as common
from mantid.api import AlgorithmFactory, DataProcessorAlgorithm, FileAction, \
    MultipleFileProperty, PropertyMode, WorkspaceGroup, WorkspaceGroupProperty
from mantid.kernel import Direction, FloatArrayProperty, FloatArrayOrderedPairsValidator, \
    FloatBoundedValidator, IntBoundedValidator, IntArrayProperty, PropertyManagerProperty, \
    RebinParamsValidator, StringListValidator
from mantid.simpleapi import *

from os import path

import numpy as np


def get_run_number(value):
    """
    Extracts the run number from the first run out of the string value of a
    multiple file property of numors
    """
    return path.splitext(path.basename(value.split(',')[0].split('+')[0]))[0]


def get_vanadium_corrections(vanadium_ws):
    """
    Extracts vanadium integral and vanadium diagnostics workspaces from the provided list. If the provided group
    has only one workspace, then it is assumed it contains vanadium integral. Assumed is the following order
    of vanadium workspaces for each numor in the group: SofQ, SofTW, diagnostics, integral.
    :param vanadium_ws: workspace group with processed vanadium
    :return: vanadium integral and vanadium diagnostics (if exist)
    """
    diagnostics = []
    integrals = []
    nentries = mtd[vanadium_ws].getNumberOfEntries()
    if nentries == 1:
        integrals.append(mtd[vanadium_ws][0].name())
    else:
        for index in range(0, nentries, 4):
            diagnostics.append(mtd[vanadium_ws][index+2].name())
            integrals.append(mtd[vanadium_ws][index+3].name())
    return diagnostics, integrals


class DirectILLAutoProcess(DataProcessorAlgorithm):

    instrument = None
    sample = None
    process = None
    reduction_type = None
    incident_energy = None
    incident_energy_calibration = None
    incident_energy_ws = None
    elastic_channel_ws = None
    masking = None
    mask_ws = None
    ebinning_params = None  # energy binning
    output = None
    vanadium = None
    vanadium_epp = None
    vanadium_diagnostics = None
    vanadium_integral = None
    empty = None
    flat_bkg_scaling = None
    flat_background = None
    to_clean = None
    absorption_corr = None
    save_output = None
    clear_cache = None

    def category(self):
        return "{};{}".format(common.CATEGORIES, "ILL\\Auto")

    def summary(self):
        return 'Performs automatic data reduction for the direct geometry TOF spectrometers at ILL.'

    def seeAlso(self):
        return ['DirectILLReduction']

    def name(self):
        return 'DirectILLAutoProcess'

    def validateInputs(self):
        issues = dict()

        run_no_err = 'Wrong number of {0} runs: {1}. Provide one or as many as '\
                     'sample runs: {2}.'
        runs_sample = len(self.getPropertyValue('Runs'))
        if not self.getProperty('EmptyContainerWorkspace').isDefault:
            runs_container = mtd[self.getPropertyValue('EmptyContainerWorkspace')].getNumberOfEntries()
            if runs_container != 1 and runs_container > runs_sample:
                issues['BeamRuns'] = run_no_err.format('EmptyContainerWorkspace', runs_container, runs_sample)

        grouping_err_msg = 'Only one grouping method can be specified.'
        if self.getProperty(common.PROP_DET_GROUPING).isDefault:
            if self.getProperty('ApplyGroupingBy').value \
                    and not self.getProperty(common.PROP_GROUPING_ANGLE_STEP).isDefault:
                issues['ApplyGroupingBy'] = grouping_err_msg
                issues[common.PROP_GROUPING_ANGLE_STEP] = grouping_err_msg
        else:
            if not self.getProperty(common.PROP_DET_GROUPING_BY).isDefault:
                issues[common.PROP_DET_GROUPING] = grouping_err_msg
                issues[common.PROP_DET_GROUPING_BY] = grouping_err_msg
            if not self.getProperty('ApplyGroupingBy').isDefault:
                issues[common.PROP_DET_GROUPING] = grouping_err_msg
                issues[common.PROP_DET_HOR_GROUPING] = grouping_err_msg
                issues[common.PROP_DET_VER_GROUPING] = grouping_err_msg
                issues['ApplyGroupingBy'] = grouping_err_msg
            if not self.getProperty(common.PROP_GROUPING_ANGLE_STEP).isDefault:
                issues[common.PROP_DET_GROUPING] = grouping_err_msg
                issues[common.PROP_GROUPING_ANGLE_STEP] = grouping_err_msg

        if not self.getProperty('VanadiumWorkspace').isDefault \
                and self.getPropertyValue('VanadiumWorkspace') not in mtd:
            # attempts to load the file, raises a runtime error if the desired file does not exist
            vanadium_ws = self.getPropertyValue('VanadiumWorkspace')
            try:
                Load(Filename=vanadium_ws,
                     OutputWorkspace=vanadium_ws)
                if not isinstance(mtd[vanadium_ws], WorkspaceGroup):
                    RenameWorkspace(InputWorkspace=vanadium_ws, OutputWorkspace='{}_1'.format(vanadium_ws))
                    GroupWorkspaces(InputWorkspaces='{}_1'.format(vanadium_ws), OutputWorkspace=vanadium_ws)
            except ValueError:
                issues['VanadiumWorkspace'] = "Desired vanadium workspace: {} cannot be found.".format(vanadium_ws)

        if self.getProperty('MaskWithVanadium').value and self.getProperty('VanadiumWorkspace').isDefault:
            issues['VanadiumWorkspace'] = 'Please provide a vanadium input for a masking reference.'

        if self.getPropertyValue('AbsorptionCorrection') != 'None':
            if self.getProperty('SampleMaterial').isDefault:
                issues['SampleMaterial'] = 'Please define sample material.'
            if self.getProperty('SampleGeometry').isDefault:
                issues['SampleGeometry'] = 'Please define sample geometry.'

        return issues

    def setUp(self):
        self.sample = self.getPropertyValue('Runs').split(',')
        self.output = self.getPropertyValue('OutputWorkspace')
        self.process = self.getPropertyValue('ProcessAs')
        self.reduction_type = self.getPropertyValue('ReductionType')
        self.to_clean = []
        if self.getProperty('IncidentEnergyCalibration').value:
            self.incident_energy_calibration = 'Energy Calibration ON'
        if not self.getProperty('IncidentEnergy').isDefault:
            self.incident_energy = self.getProperty('IncidentEnergy').value
            self.incident_energy_ws = 'incident_energy_ws'
            CreateSingleValuedWorkspace(DataValue=self.incident_energy,
                                        OutputWorkspace=self.incident_energy_ws)
            self.to_clean.append(self.incident_energy_ws)
        if not self.getProperty('ElasticChannelIndex').isDefault:
            self.elastic_channel_ws = 'elastic_channel_index_ws'
            CreateSingleValuedWorkspace(DataValue=self.getProperty('ElasticChannelIndex').value,
                                        OutputWorkspace=self.elastic_channel_ws)
            self.to_clean.append(self.elastic_channel_ws)
        if (self.getProperty('MaskWorkspace').isDefault and self.getProperty('MaskedTubes').isDefault
                and self.getProperty('MaskThresholdMin').isDefault and self.getProperty('MaskThresholdMax').isDefault
                and self.getProperty('MaskedAngles').isDefault and self.getProperty('MaskWithVanadium').isDefault):
            self.masking = False
        else:
            self.masking = True
        self.flat_bkg_scaling = self.getProperty(common.PROP_FLAT_BKG_SCALING).value
        self.ebinning_params = self.getProperty('EnergyExchangeBinning').value
        self.empty = self.getPropertyValue('EmptyContainerWorkspace')
        self.vanadium = self.getPropertyValue('VanadiumWorkspace')
        if self.vanadium != str():
            self.vanadium_diagnostics, self.vanadium_integral = get_vanadium_corrections(self.vanadium)
        self.save_output = self.getProperty('SaveOutput').value
        self.clear_cache = self.getProperty('ClearCache').value

    def PyInit(self):

        positiveFloat = FloatBoundedValidator(0., exclusive=False)
        positiveInt = IntBoundedValidator(0, exclusive=False)
        validRebinParams = RebinParamsValidator(AllowEmpty=True)
        orderedPairsValidator = FloatArrayOrderedPairsValidator()

        self.declareProperty(WorkspaceGroupProperty('OutputWorkspace', '',
                                                    direction=Direction.Output),
                             doc='The output workspace group containing reduced data.')

        self.declareProperty(MultipleFileProperty('Runs',
                                                  action=FileAction.Load,
                                                  extensions=['nxs']),
                             doc='Run(s) to be processed.')

        processes = ['Cadmium', 'Empty', 'Vanadium', 'Sample']
        self.declareProperty(name='ProcessAs',
                             defaultValue='Sample',
                             validator=StringListValidator(processes),
                             doc='Choose the process type.')

        reduction_options = ['Powder', 'SingleCrystal']
        self.declareProperty(name='ReductionType',
                             defaultValue='Powder',
                             validator=StringListValidator(reduction_options),
                             doc='Choose the appropriate reduction type for the data to process.')

        self.declareProperty('VanadiumWorkspace', '',
                             doc='File(s) or workspaces containing vanadium data.')

        self.declareProperty(WorkspaceGroupProperty('EmptyContainerWorkspace', '',
                                                    direction=Direction.Input,
                                                    optional=PropertyMode.Optional),
                             doc='Empty container workspace.')

        self.declareProperty('EmptyContainerScaling', 1.0,
                             doc='Scaling factor for the empty container.')

        self.declareProperty(WorkspaceGroupProperty('CadmiumWorkspace', '',
                                                    direction=Direction.Input,
                                                    optional=PropertyMode.Optional),
                             doc='Cadmium absorber workspace.')

        self.copyProperties('DirectILLCollectData', [common.PROP_FLAT_BKG, common.PROP_FLAT_BKG_WINDOW])

        self.declareProperty('FlatBackgroundSource', "",
                             doc='File(s) or workspaces containing the source to calculate flat background.')

        self.copyProperties('DirectILLCollectData', common.PROP_FLAT_BKG_SCALING)

        self.declareProperty(common.PROP_ABSOLUTE_UNITS, False,
                             doc='Enable or disable normalisation to absolute units.')

        self.declareProperty("IncidentEnergyCalibration", True,
                             doc='Enable or disable incident energy calibration.')

        self.declareProperty("ElasticChannelCalibration", False,
                             doc='Enable or disable calibration of the elastic channel.')

        additional_inputs_group = 'Corrections'
        self.setPropertyGroup('VanadiumWorkspace', additional_inputs_group)
        self.setPropertyGroup('EmptyContainerWorkspace', additional_inputs_group)
        self.setPropertyGroup('EmptyContainerScaling', additional_inputs_group)
        self.setPropertyGroup('CadmiumWorkspace', additional_inputs_group)
        self.setPropertyGroup(common.PROP_FLAT_BKG, additional_inputs_group)
        self.setPropertyGroup(common.PROP_FLAT_BKG_WINDOW, additional_inputs_group)
        self.setPropertyGroup('FlatBackgroundSource', additional_inputs_group)
        self.setPropertyGroup(common.PROP_FLAT_BKG, additional_inputs_group)
        self.setPropertyGroup('IncidentEnergyCalibration', additional_inputs_group)
        self.setPropertyGroup('ElasticChannelCalibration', additional_inputs_group)
        self.setPropertyGroup(common.PROP_ABSOLUTE_UNITS, additional_inputs_group)

        self.declareProperty(name='IncidentEnergy',
                             defaultValue=0.0,
                             validator=positiveFloat,
                             doc='Value for the calibrated incident energy (meV).')

        self.declareProperty(name='ElasticChannelIndex',
                             defaultValue=0,
                             validator=positiveInt,
                             doc='Index number for the centre of the elastic peak.')

        self.declareProperty('SampleAngleOffset', 0.0,
                             doc='Value for the offset parameter in omega scan (degrees).')

        parameters_group = 'Parameters'
        self.setPropertyGroup('IncidentEnergy', parameters_group)
        self.setPropertyGroup('ElasticChannelIndex', parameters_group)
        self.setPropertyGroup('SampleAngleOffset', parameters_group)

        # The mask workspace replaces MaskWorkspace parameter from PantherSingle and DiagnosticsWorkspace from directred
        self.declareProperty('MaskWorkspace', '',
                             doc='File(s) or workspaces containing the mask.')

        self.declareProperty(IntArrayProperty(name='MaskedTubes', values=[], direction=Direction.Input),
                             doc='List of tubes to be masked.')

        self.declareProperty('MaskThresholdMin', 0.0, doc='Threshold level below which bins will be masked'
                                                          ' to remove empty / background pixels.')

        self.declareProperty('MaskThresholdMax', 0.0, doc='Threshold level above which bins will be masked'
                                                          ' to remove noisy pixels.')

        self.declareProperty(FloatArrayProperty(name='MaskedAngles', values=[], validator=orderedPairsValidator),
                             doc='Mask detectors in the given angular range.')

        self.declareProperty('MaskWithVanadium', False,
                             doc='Whether to mask using vanadium workspace.')

        masking_group_name = 'Masking'
        self.setPropertyGroup('MaskWorkspace', masking_group_name)
        self.setPropertyGroup('MaskedTubes', masking_group_name)
        self.setPropertyGroup('MaskThresholdMin', masking_group_name)
        self.setPropertyGroup('MaskThresholdMax', masking_group_name)
        self.setPropertyGroup('MaskedAngles', masking_group_name)
        self.setPropertyGroup('MaskWithVanadium', masking_group_name)

        self.declareProperty(FloatArrayProperty(name='EnergyExchangeBinning',
                                                validator=validRebinParams),
                             doc='Energy exchange binning parameters.')

        self.declareProperty(FloatArrayProperty(name='MomentumTransferBinning',
                                                validator=validRebinParams),
                             doc='Momentum transfer binning parameters.')

        rebinning_group = 'Binning parameters'
        self.setPropertyGroup('EnergyExchangeBinning', rebinning_group)
        self.setPropertyGroup('MomentumTransferBinning', rebinning_group)

        self.declareProperty(name='AbsorptionCorrection',
                             defaultValue='None',
                             validator=StringListValidator(['None', 'Fast', 'Full']),
                             doc='Choice of approach to absorption correction.')

        self.declareProperty(name='SelfAttenuationMethod',
                             defaultValue='MonteCarlo',
                             validator=StringListValidator(['Numerical', 'MonteCarlo']),
                             doc='Choice of calculation method for the attenuation calculation.')

        self.declareProperty(PropertyManagerProperty('SampleMaterial'),
                             doc='Sample material definitions.')

        self.declareProperty(PropertyManagerProperty('SampleGeometry'),
                             doc="Dictionary for the sample geometry.")

        self.declareProperty(PropertyManagerProperty('ContainerMaterial'),
                             doc='Container material definitions.')

        self.declareProperty(PropertyManagerProperty('ContainerGeometry'),
                             doc="Dictionary for the container geometry.")

        attenuation_group = 'Sample attenuation'
        self.setPropertyGroup('AbsorptionCorrection', attenuation_group)
        self.setPropertyGroup('SelfAttenuationMethod', attenuation_group)
        self.setPropertyGroup('SampleMaterial', attenuation_group)
        self.setPropertyGroup('SampleGeometry', attenuation_group)
        self.setPropertyGroup('ContainerMaterial', attenuation_group)
        self.setPropertyGroup('ContainerGeometry', attenuation_group)

        self.declareProperty(name=common.PROP_DET_GROUPING,
                             defaultValue="",
                             doc='Grouping pattern to reduce the granularity of the output.')

        self.declareProperty(name=common.PROP_DET_GROUPING_BY,
                             defaultValue=1,
                             doc='Step to use when grouping detectors to reduce the granularity of the output.')

        self.copyProperties('DirectILLCollectData', [common.PROP_DET_HOR_GROUPING, common.PROP_DET_VER_GROUPING])

        self.declareProperty(name="ApplyGroupingBy",
                             defaultValue=False,
                             doc='Whether to apply the pixel grouping horizontally or vertically to the data, and not'
                                 ' only to increase the statistics of the flat background calculation.')

        self.declareProperty(name=common.PROP_GROUPING_ANGLE_STEP,
                             defaultValue=0.0,
                             validator=positiveFloat,
                             doc='A scattering angle step to which to group detectors, in degrees.')

        self.declareProperty(name='GroupingBehaviour',
                             defaultValue="Sum",
                             validator=StringListValidator(['Sum', 'Average']),
                             doc='Defines which behaviour should be used when grouping pixels.')

        grouping_options_group = 'Grouping options'
        self.setPropertyGroup(common.PROP_DET_GROUPING, grouping_options_group)
        self.setPropertyGroup(common.PROP_DET_GROUPING_BY, grouping_options_group)
        self.setPropertyGroup(common.PROP_DET_HOR_GROUPING, grouping_options_group)
        self.setPropertyGroup(common.PROP_DET_VER_GROUPING, grouping_options_group)
        self.setPropertyGroup('ApplyGroupingBy', grouping_options_group)
        self.setPropertyGroup(common.PROP_GROUPING_ANGLE_STEP, grouping_options_group)
        self.setPropertyGroup('GroupingBehaviour', grouping_options_group)

        self.declareProperty(name="SaveOutput",
                             defaultValue=True,
                             doc="Whether to save the output directly after processing.")

        self.declareProperty(name='ClearCache',
                             defaultValue=False,
                             doc='Whether to clear intermediate workspaces.')

    def PyExec(self):
        self.setUp()
        sample_runs = self.getPropertyValue('Runs').split(',')
        output_samples = []
        for sample_no, sample in enumerate(sample_runs):
            current_output = []  # output of the current iteration of reduction
            ws = self._collect_data(sample, vanadium=self.process == 'Vanadium')
            if self.masking and sample_no == 0:  # prepares masks once, and when the instrument is known
                self.mask_ws = self._prepare_masks()
            if self.process == 'Vanadium':
                ws_sofq, ws_softw, ws_diag, ws_integral = self._process_vanadium(ws)
                current_output = [ws_sofq, ws_softw, ws_diag, ws_integral]
                output_samples.extend(current_output)
            elif self.process == 'Sample':
                sample_sofq, sample_softw = self._process_sample(ws, sample_no)
                current_output = np.array([sample_sofq, sample_softw])
                current_output = current_output[[isinstance(elem, str) for elem in current_output]]
                output_samples.extend(current_output)
            else:  # Empty or Cadmium
                current_output = ws
                output_samples.append(current_output)
            self._group_detectors(current_output)
            if self.save_output:
                self._save_output(current_output)
        GroupWorkspaces(InputWorkspaces=output_samples,
                        OutputWorkspace=self.output)
        if self.clear_cache:  # final clean up
            self._clean_up(self.to_clean)
        self.setProperty('OutputWorkspace', mtd[self.output])

    def _collect_data(self, sample, vanadium=False):
        """Loads data if the corresponding workspace does not exist in the ADS."""
        ws = "{}_{}".format(get_run_number(sample), 'raw')
        kwargs = dict()
        if not self.getProperty('FlatBackgroundSource').isDefault:
            kwargs['FlatBkgWorkspace'] = self.getPropertyValue('FlatBackgroundSource')
            kwargs[common.PROP_FLAT_BKG_SCALING] = self.flat_bkg_scaling
        elif not self.getProperty(common.PROP_FLAT_BKG).isDefault:
            kwargs[common.PROP_FLAT_BKG] = self.getPropertyValue(common.PROP_FLAT_BKG)
            kwargs[common.PROP_FLAT_BKG_WINDOW] = self.getPropertyValue(common.PROP_FLAT_BKG_WINDOW)
        if not self.getProperty(common.PROP_DET_HOR_GROUPING).isDefault:
            kwargs[common.PROP_DET_HOR_GROUPING] = self.getProperty(common.PROP_DET_HOR_GROUPING).value
        if not self.getProperty(common.PROP_DET_VER_GROUPING).isDefault:
            kwargs[common.PROP_DET_VER_GROUPING] = self.getProperty(common.PROP_DET_VER_GROUPING).value
        if vanadium:
            kwargs['EPPCreationMethod'] = 'Calculate EPP'
            kwargs['ElasticChannel'] = 'Elastic Channel AUTO'
            kwargs['FlatBkg'] = 'Flat Bkg ON'
            self.vanadium_epp = "{}_epp".format(ws)
            kwargs['OutputEPPWorkspace'] = self.vanadium_epp
            self.to_clean.append(self.vanadium_epp)
        DirectILLCollectData(Run=sample, OutputWorkspace=ws,
                             IncidentEnergyCalibration=self.incident_energy_calibration,
                             IncidentEnergyWorkspace=self.incident_energy_ws,
                             ElasticChannelWorkspace=self.elastic_channel_ws,
                             **kwargs)
        instrument = mtd[ws].getInstrument().getName()
        if self.instrument and instrument != self.instrument:
            self.log().error("Sample data: {} comes from different instruments that the rest of the data:"
                             " {} and {}".format(sample, instrument, self.instrument))
        else:
            self.instrument = instrument
        return ws

    @staticmethod
    def _clean_up(to_clean):
        """Performs the clean up of intermediate workspaces that are created and used throughout the code."""
        if len(to_clean) > 0:
            DeleteWorkspaces(WorkspaceList=to_clean)

    def _save_output(self, ws_to_save):
        """Saves the output workspaces to an external file."""
        for ws_name in ws_to_save:
            if self.reduction_type == 'SingleCrystal':
                omega_log = 'SRot.value'  # IN5, IN6
                if self.instrument == 'PANTHER':
                    omega_log = 'a3.value'
                elif self.instrument == 'SHARP':
                    omega_log = 'srotation.value'
                elif self.instrument == 'IN4':
                    omega_log = 'SampleRotation.value'
                psi = mtd[ws_name].run().getProperty(omega_log).value
                psi_offset = self.getProperty('SampleAngleOffset').value
                offset = psi - psi_offset
                SaveNXSPE(
                    InputWorkspace=ws_name,
                    Filename='{}.nxspe'.format(ws_name),
                    Psi=offset
                )
            else:  # powder reduction
                SaveNexus(
                    InputWorkspace=ws_name,
                    Filename='{}.nxs'.format(ws_name)
                )

    def _get_grouping_pattern(self, ws):
        """Returns a grouping pattern taking into account the requested grouping, either inside the tube, as defined
         by grouping-by property, or horizontally (between tubes) and vertically (inside a tube).
         The latter method can be applied only to PANTHER, SHARP, and IN5 instruments' structure."""
        n_pixels = mtd[ws].getNumberHistograms()
        if not self.getProperty(common.PROP_DET_GROUPING_BY).isDefault:
            group_by = self.getProperty(common.PROP_DET_GROUPING_BY).value
            grouping_pattern = \
                ["{}-{}".format(pixel_id, pixel_id + group_by - 1) for pixel_id in range(0, n_pixels, group_by)]
        else:
            group_by_x = self.getProperty(common.PROP_DET_HOR_GROUPING).value
            group_by_y = self.getProperty(common.PROP_DET_VER_GROUPING).value
            n_tubes = 0
            n_monitors = {'IN5': 1, 'PANTHER': 1, 'SHARP': 1}
            for comp in mtd[ws].componentInfo():
                if len(comp.detectorsInSubtree) > 1 and comp.hasParent:
                    n_tubes += 1
            n_tubes -= n_monitors[self.instrument]
            if self.instrument == 'IN5':  # there is an extra bank that contains all of the tubes
                n_tubes -= 1
            n_pixels_per_tube = int(n_pixels / n_tubes)
            print("N tubes:", n_tubes, 'npix/tube:', n_pixels_per_tube)
            grouping_pattern = []
            pixel_id = 0
            while pixel_id < n_pixels - (group_by_x - 1) * n_pixels_per_tube:
                pattern = []
                for tube_shift in range(0, group_by_x):
                    numeric_pattern = list(range(pixel_id + tube_shift * n_pixels_per_tube,
                                                 pixel_id + tube_shift * n_pixels_per_tube + group_by_y))
                    pattern.append('+'.join(map(str, numeric_pattern)))
                pattern = "+".join(pattern)
                grouping_pattern.append(pattern)
                pixel_id += group_by_y
                if pixel_id % n_pixels_per_tube == 0:
                    pixel_id += n_pixels_per_tube * (group_by_x - 1)

        return ",".join(grouping_pattern)

    def _group_detectors(self, ws_list):
        """Groups detectors for workspaces in the provided list according to the defined grouping pattern."""
        grouping_pattern = None
        if not self.getProperty(common.PROP_DET_GROUPING).isDefault:
            grouping_pattern = self.getProperty(common.PROP_DET_GROUPING).value
        elif not self.getProperty(common.PROP_DET_GROUPING_BY).isDefault or self.getProperty('ApplyGroupingBy').value:
            grouping_pattern = self._get_grouping_pattern(ws_list[0])
        if grouping_pattern is not None:
            for ws in ws_list:
                GroupDetectors(InputWorkspace=ws, OutputWorkspace=ws,
                               GroupingPattern=grouping_pattern,
                               Behaviour=self.getPropertyValue('GroupingBehaviour'))

    def _prepare_masks(self):
        """Builds a masking workspace from the provided inputs. Masking using threshold cannot be prepared ahead."""
        existing_masks = []
        mask = self.getPropertyValue('MaskWorkspace')
        if mask != str():
            mask = self.getPropertyValue('MaskWorkspace')
            if mask not in mtd:
                LoadMask(Instrument=self.instrument, InputFile=mask, OutputWorkspace=mask)
            existing_masks.append(mask)
        mask_tubes = self.getPropertyValue('MaskedTubes')
        if mask_tubes != str():
            MaskBTP(Instrument=self.instrument, Tube=self.getPropertyValue(mask_tubes))
            tube_mask_ws = "{}_masked_tubes".format(self.instrument)
            RenameWorkspace(InputWorkspace='{}MaskBTP'.format(self.instrument), OutputWorkspace=tube_mask_ws)
            existing_masks.append(tube_mask_ws)

        mask_angles = self.getProperty('MaskedAngles').value
        if mask_angles != list():
            masked_angles_ws = '{}_masked_angles'.format(self.instrument)
            LoadEmptyInstrument(Filename=self.instrument, OutputWorkspace=masked_angles_ws)
            MaskAngles(Workspace=masked_angles_ws, MinAngle=mask_angles[0], MaxAngle=mask_angles[1])
            existing_masks.append(masked_angles_ws)
        mask_with_vanadium = self.getProperty('MaskWithVanadium').value
        if mask_with_vanadium:
            existing_masks.append(self.vanadium_diagnostics)

        mask_ws = 'mask_ws'
        if len(existing_masks) > 1:
            MergeRuns(InputWorkspaces=existing_masks, OutputWorkspace=mask_ws)
        else:
            RenameWorkspace(InputWorkspace=existing_masks[0], OutputWorkspace=mask_ws)
        self.to_clean.append(mask_ws)
        return mask_ws

    def _apply_mask(self, ws):
        """Applies selected masks."""
        MaskDetectors(Workspace=ws, MaskedWorkspace=self.mask_ws)
        # masks bins below the chosen threshold, this has to be applied for each ws and cannot be created ahead:
        min_threshold_defined = not self.getProperty('MaskThresholdMin').isDefault
        max_threshold_defined = not self.getProperty('MaskThresholdMax').isDefault
        if not min_threshold_defined or not max_threshold_defined:
            masking_criterion = '{} < y < {}'
            if min_threshold_defined and max_threshold_defined:
                masking_criterion = masking_criterion.format(self.getPropertyValue('MaskThresholdMax'),
                                                             self.getPropertyValue('MaskThresholdMin'))
            elif min_threshold_defined:
                masking_criterion = 'y < {}'.format(self.getPropertyValue('MaskThresholdMin'))
            else:  # only max threshold defined
                masking_criterion = 'y > {}'.format(self.getPropertyValue('MaskThresholdMax'))

            MaskBinsIf(InputWorkspace=ws, OutputWorkspace=ws,
                       Criterion=masking_criterion)
        return ws

    def _subtract_empty_container(self, ws):
        """Subtracts empty container counts from the sample workspace."""
        empty_ws = self.getPropertyValue('EmptyContainerWorkspace')
        empty_correction_ws = "{}_correction".format(empty_ws)
        empty_scaling = self.getProperty('EmptyContainerScaling').value
        Scale(InputWorkspace=empty_ws,
              OutputWorkspace=empty_correction_ws,
              Factor=empty_scaling)
        Minus(LHSWorkspace=ws,
              RHSWorkspace=mtd[empty_correction_ws][0],
              OutputWorkspace=ws)
        if self.clear_cache:
            DeleteWorkspace(Workspace=empty_correction_ws)

    def _process_vanadium(self, ws):
        """Processes vanadium and creates workspaces with diagnostics, integrated vanadium, and reduced vanadium."""
        to_remove = [ws]
        numor = ws[:ws.rfind('_')]
        vanadium_diagnostics = '{}_diag'.format(numor)
        DirectILLDiagnostics(InputWorkspace=ws,
                             OutputWorkspace=vanadium_diagnostics,
                             BeamStopDiagnostics="Beam Stop Diagnostics OFF")

        if self.instrument == 'IN5':
            van_flat_ws = "{}_flat".format(numor)
            to_remove.append(van_flat_ws)
            DirectILLTubeBackground(InputWorkspace=ws,
                                    DiagnosticsWorkspace=vanadium_diagnostics,
                                    EPPWorkspace=self.vanadium_epp,
                                    OutputWorkspace=van_flat_ws)
            Minus(LHSWorkspace=ws, RHSWorkspace=van_flat_ws,
                  OutputWorkspace=ws, EnableLogging=False)

        if self.empty:
            self._subtract_empty_container(ws)

        vanadium_integral = '{}_integral'.format(numor)
        DirectILLIntegrateVanadium(InputWorkspace=ws,
                                   OutputWorkspace=vanadium_integral,
                                   EPPWorkspace=self.vanadium_epp)

        sofq_output = '{}_SofQ'.format(numor)
        softw_output = '{}_SofTW'.format(numor)
        optional_parameters = dict()
        if not self.getProperty(common.PROP_GROUPING_ANGLE_STEP).isDefault:
            optional_parameters['GroupingAngleStep'] = self.getProperty(common.PROP_GROUPING_ANGLE_STEP).value
        if not self.getProperty('EnergyExchangeBinning').isDefault:
            optional_parameters['EnergyRebinningParams'] = self.getProperty(common.PROP_GROUPING_ANGLE_STEP).value
        if not self.getProperty('MomentumTransferBinning').isDefault:
            optional_parameters['QBinningParams'] = self.getProperty(MomentumTransferBinning).value
        DirectILLReduction(InputWorkspace=ws,
                           OutputWorkspace=sofq_output,
                           OutputSofThetaEnergyWorkspace=softw_output,
                           IntegratedVanadiumWorkspace=vanadium_integral,
                           DiagnosticsWorkspace=vanadium_diagnostics,
                           **optional_parameters
                           )

        if len(to_remove) > 0 and self.clear_cache:
            self._clean_up(to_remove)
        return sofq_output, softw_output, vanadium_diagnostics, vanadium_integral

    def _normalise_sample(self, sample_ws, sample_no, numor):
        """
        Normalises sample using vanadium integral, if it has been provided.
        :param sample_ws: sample being processed
        :return: Either normalised sample or the input, if vanadium is not provided
        """
        normalised_ws = '{}_norm'.format(numor)
        if self.vanadium_integral and self.vanadium_integral != list():
            nintegrals = len(self.vanadium_integral)
            vanadium_no = sample_no
            if nintegrals == 1:
                vanadium_no = 0
            elif sample_no > nintegrals:
                vanadium_no = sample_no % nintegrals
            Divide(
                LHSWorkspace=sample_ws,
                RHSWorkspace=self.vanadium_integral[vanadium_no],
                OutputWorkspace=normalised_ws
            )
        else:
            normalised_ws = sample_ws
            self.log().warning("Vanadium integral workspace not found.\nData will not be normalised to vanadium.")

        return normalised_ws

    def _prepare_self_attenuation_ws(self, ws):
        """Creates a self-attenuation workspace using either a MonteCarlo approach or numerical integration."""
        sample_geometry = self.getProperty('SampleGeometry').value
        sample_material = self.getProperty('SampleMaterial').value
        container_geometry = self.getProperty('ContainerGeometry').value \
            if not self.getProperty('SampleGeometry').isDefault else ""
        container_material = self.getProperty('ContainerMaterial').value \
            if not self.getProperty('ContainerMaterial').isDefault else ""
        self.absorption_corr = "{}_abs_corr".format(ws)
        self.to_clean.append(self.absorption_corr)
        SetSample(
            InputWorkspace=ws,
            Geometry=sample_geometry,
            Material=sample_material,
            ContainerGeometry=container_geometry,
            ContainerMaterial=container_material
        )
        if self.getProperty('SelfAttenuationMethod').value == 'MonteCarlo':
            sparse_parameters = dict()
            if self.getPropertyValue('AbsorptionCorrection') == 'Fast':
                sparse_parameters['SparseInstrument'] = True
                n_detector_rows = 5
                n_detector_cols = 10
                sparse_parameters['NumberOfDetectorRows'] = n_detector_rows
                sparse_parameters['NumberOfDetectorColumns'] = n_detector_cols

            PaalmanPingsMonteCarloAbsorption(
                InputWorkspace=ws,
                CorrectionsWorkspace=self.absorption_corr,
                **sparse_parameters
            )
        else:
            PaalmanPingsAbsorptionCorrection(
                InputWorkspace=ws,
                OutputWorkspace=self.absorption_corr
            )

    def _correct_self_attenuation(self, ws, sample_no):
        """Creates, if necessary, a self-attenuation workspace and uses it to correct the provided sample workspace."""
        if sample_no == 0:
            self._prepare_self_attenuation_ws(ws)

        if self.absorption_corr:
            ApplyPaalmanPingsCorrection(
                    SampleWorkspace=ws,
                    OutputWorkspace=ws,
                    CorrectionsWorkspace=self.absorption_corr
            )

    def _process_sample(self, ws, sample_no):
        """Does the sample data reduction for single crystal."""
        to_remove = [ws]
        if self.empty:
            self._subtract_empty_container(ws)
        if self.masking:
            ws = self._apply_mask(ws)
        numor = ws[:ws.rfind('_')]
        processed_sample_tw = None
        if self.reduction_type == 'SingleCrystal':
            # normalises to vanadium integral
            normalised_ws = self._normalise_sample(ws, sample_no, numor)
            to_remove.append(normalised_ws)
            # converts to energy
            corrected_ws = '{}_ene'.format(numor)
            ConvertUnits(InputWorkspace=normalised_ws, EFixed=self.incident_energy,
                         Target='DeltaE', EMode='Direct', OutputWorkspace=corrected_ws)
            to_remove.append(corrected_ws)

            # transforms the distribution into dynamic structure factor
            CorrectKiKf(InputWorkspace=corrected_ws, EFixed=self.incident_energy,
                        OutputWorkspace=corrected_ws)

            # corrects for detector efficiency
            DetectorEfficiencyCorUser(InputWorkspace=corrected_ws, IncidentEnergy=self.incident_energy,
                                      OutputWorkspace=corrected_ws)

            # rebin in energy or momentum transfer
            processed_sample = '{}_reb'.format(numor)
            if self.ebinning_params != list():
                Rebin(InputWorkspace=corrected_ws, Params=self.ebinning_params, OutputWorkspace=processed_sample)
            else:
                RenameWorkspace(InputWorkspace=corrected_ws, OutputWorkspace=processed_sample)
                to_remove.pop()
            # saving of the output is omitted at this point: it is handled by Drill interface
        else:
            processed_sample = 'SofQW_{}'.format(ws[:ws.rfind('_')])  # name should contain only SofQW and numor
            processed_sample_tw = 'SofTW_{}'.format(ws[:ws.rfind('_')])  # name should contain only SofTW and numor
            if self.getPropertyValue('AbsorptionCorrection') != 'None':
                self._correct_self_attenuation(ws, sample_no)
            vanadium_integral = self.vanadium_integral[0] if self.vanadium_integral else ""
            vanadium_diagnostics = self.vanadium_diagnostics[0] if self.vanadium_diagnostics else ""
            optional_parameters = dict()
            if not self.getProperty(common.PROP_GROUPING_ANGLE_STEP).isDefault:
                optional_parameters['GroupingAngleStep'] = self.getProperty(common.PROP_GROUPING_ANGLE_STEP).value
            if not self.getProperty('EnergyExchangeBinning').isDefault:
                optional_parameters['EnergyRebinningParams'] = self.getProperty(common.PROP_GROUPING_ANGLE_STEP).value
            if not self.getProperty('MomentumTransferBinning').isDefault:
                optional_parameters['QBinningParams'] = self.getProperty(MomentumTransferBinning).value
            DirectILLReduction(
                InputWorkspace=ws,
                OutputWorkspace=processed_sample,
                OutputSofThetaEnergyWorkspace=processed_sample_tw,
                IntegratedVanadiumWorkspace=vanadium_integral,
                DiagnosticsWorkspace=vanadium_diagnostics,
                **optional_parameters
            )
        if len(to_remove) > 0 and self.clear_cache:
            self._clean_up(to_remove)

        return processed_sample, processed_sample_tw


AlgorithmFactory.subscribe(DirectILLAutoProcess)
