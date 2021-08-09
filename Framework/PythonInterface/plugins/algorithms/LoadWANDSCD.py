# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
from mantid.api import (PythonAlgorithm, AlgorithmFactory, PropertyMode, ReplicateMD, FileProperty, DivideMD, WorkspaceProperty, Progress,
                        MultipleFileProperty, FileAction, mtd)
from mantid.kernel import Direction, Property, IntArrayProperty, StringListValidator, FloatTimeSeriesProperty
from mantid.simpleapi import LoadEventNexus, RemoveLogs, DeleteWorkspace, ConvertToMD, Rebin, CreateGroupingWorkspace, GroupDetectors, SetUB
import numpy as np
import h5py
import re


class LoadWANDSCD(PythonAlgorithm):
    def category(self):
        return 'DataHandling\\Nexus'

    def seeAlso(self):
        return ["ConvertWANDSCDtoQ"]

    def name(self):
        return 'LoadWANDSCD'

    def summary(self):
        return 'Load WAND single crystal data into a detector space vs rotation MDHisto'

    def PyInit(self):
        # Input workspace/data info (filename or IPTS+RunNumber)
        self.declareProperty(MultipleFileProperty(name="Filename", action=FileAction.OptionalLoad, extensions=[".nxs.h5"]), "Files to load")
        self.declareProperty('IPTS', Property.EMPTY_INT, "IPTS number to load from")
        self.declareProperty(IntArrayProperty("RunNumbers", []), 'Run numbers to load')
        # Normalization info (optional, skip normalization if not specified)
        self.declareProperty(FileProperty(name="NormWorkspace", action=FileAction.OptionalLoad, extensions=[".nxs.h5"]),
                             "Normalization workspace")
        self.declareProperty("NormalizedBy", '', StringListValidator(['', 'counts', 'Monitor', 'Time']), "Normalisation")
        # Output workspace/data info
        self.declareProperty("Grouping", 'None', StringListValidator(['None', '2x2', '4x4']), "Group pixels")
        self.declareProperty(WorkspaceProperty("OutputWorkspace", "", optional=PropertyMode.Mandatory, direction=Direction.Output),
                             "Output Workspace")

    def validateInputs(self):
        issues = dict()

        if not self.getProperty("Filename").value:
            if (self.getProperty("IPTS").value == Property.EMPTY_INT) or len(self.getProperty("RunNumbers").value) == 0:
                issues["Filename"] = 'Must specify either Filename or IPTS AND RunNumbers'

        return issues

    def PyExec(self):
        runs = self.getProperty("Filename").value

        if not runs:
            ipts = self.getProperty("IPTS").value
            runs = ['/HFIR/HB2C/IPTS-{}/nexus/HB2C_{}.nxs.h5'.format(ipts, run) for run in self.getProperty("RunNumbers").value]

        grouping = self.getProperty("Grouping").value
        if grouping == 'None':
            grouping = 1
        else:
            grouping = 2 if grouping == '2x2' else 4

        x_dim = 480*8 // grouping
        y_dim = 512 // grouping

        number_of_runs = len(runs)

        data_array = np.empty((number_of_runs, x_dim, y_dim), dtype=np.float64)

        s1_array = []
        duration_array = []
        run_number_array = []
        monitor_count_array = []

        progress = Progress(self, 0.0, 1.0, number_of_runs+3)

        for n, run in enumerate(runs):
            progress.report('Loading: '+run)
            with h5py.File(run, 'r') as f:
                bc = np.zeros((512*480*8),dtype=np.int64)
                for b in range(8):
                    bc += np.bincount(f['/entry/bank'+str(b+1)+'_events/event_id'].value,minlength=512*480*8)
                bc = bc.reshape((480*8, 512))
                if grouping == 2:
                    bc = bc[::2,::2]+bc[1::2,::2]+bc[::2,1::2]+bc[1::2,1::2]
                elif grouping == 4:
                    bc = (bc[::4,::4]    + bc[1::4,::4]  + bc[2::4,::4]  + bc[3::4,::4]
                          + bc[::4,1::4] + bc[1::4,1::4] + bc[2::4,1::4] + bc[3::4,1::4]
                          + bc[::4,2::4] + bc[1::4,2::4] + bc[2::4,2::4] + bc[3::4,2::4]
                          + bc[::4,3::4] + bc[1::4,3::4] + bc[2::4,3::4] + bc[3::4,3::4])
                data_array[n] = bc
                s1_array.append(f['/entry/DASlogs/HB2C:Mot:s1.RBV/average_value'].value[0])
                duration_array.append(float(f['/entry/duration'].value[0]))
                run_number_array.append(float(f['/entry/run_number'].value[0]))
                monitor_count_array.append(float(f['/entry/monitor1/total_counts'].value[0]))

        progress.report('Creating MDHistoWorkspace')
        createWS_alg = self.createChildAlgorithm("CreateMDHistoWorkspace", enableLogging=False)
        createWS_alg.setProperty("SignalInput", data_array)
        createWS_alg.setProperty("ErrorInput", np.sqrt(data_array))
        createWS_alg.setProperty("Dimensionality", 3)
        createWS_alg.setProperty("Extents", '0.5,{},0.5,{},0.5,{}'.format(y_dim+0.5, x_dim+0.5, number_of_runs+0.5))
        createWS_alg.setProperty("NumberOfBins", '{},{},{}'.format(y_dim,x_dim,number_of_runs))
        createWS_alg.setProperty("Names", 'y,x,scanIndex')
        createWS_alg.setProperty("Units", 'bin,bin,number')
        createWS_alg.execute()
        outWS = createWS_alg.getProperty("OutputWorkspace").value

        progress.report('Getting IDF')
        # Get the instrument and some logs from the first file; assume the rest are the same
        _tmp_ws = LoadEventNexus(runs[0], MetaDataOnly=True, EnableLogging=False)
        # The following logs should be the same for all runs
        RemoveLogs(_tmp_ws,
                   KeepLogs='HB2C:Mot:detz,HB2C:Mot:detz.RBV,HB2C:Mot:s2,HB2C:Mot:s2.RBV,'
                   'HB2C:Mot:sgl,HB2C:Mot:sgl.RBV,HB2C:Mot:sgu,HB2C:Mot:sgu.RBV,'
                   'run_title,start_time,experiment_identifier,HB2C:CS:CrystalAlign:UBMatrix',
                   EnableLogging=False)

        time_ns_array = _tmp_ws.run().startTime().totalNanoseconds() + np.append(0, np.cumsum(duration_array)*1e9)[:-1]

        try:
            ub = np.array(re.findall(r'-?\d+\.*\d*', _tmp_ws.run().getProperty('HB2C:CS:CrystalAlign:UBMatrix').value[0]),
                          dtype=np.float).reshape(3,3)
            sgl = np.deg2rad(_tmp_ws.run().getProperty('HB2C:Mot:sgl.RBV').value[0]) # 'HB2C:Mot:sgl.RBV,1,0,0,-1'
            sgu = np.deg2rad(_tmp_ws.run().getProperty('HB2C:Mot:sgu.RBV').value[0]) # 'HB2C:Mot:sgu.RBV,0,0,1,-1'
            sgl_a = np.array([[           1,            0,           0],
                              [           0,  np.cos(sgl), np.sin(sgl)],
                              [           0, -np.sin(sgl), np.cos(sgl)]])
            sgu_a = np.array([[ np.cos(sgu),  np.sin(sgu),           0],
                              [-np.sin(sgu),  np.cos(sgu),           0],
                              [           0,            0,           1]])
            UB = sgl_a.dot(sgu_a).dot(ub) # Apply the Goniometer tilts to the UB matrix
            SetUB(_tmp_ws, UB=UB, EnableLogging=False)
        except (RuntimeError, ValueError):
            SetUB(_tmp_ws, EnableLogging=False)

        if grouping > 1:
            _tmp_group, _, _ = CreateGroupingWorkspace(InputWorkspace=_tmp_ws, EnableLogging=False)

            group_number = 0
            for x in range(0,480*8,grouping):
                for y in range(0,512,grouping):
                    group_number += 1
                    for j in range(grouping):
                        for i in range(grouping):
                            _tmp_group.dataY(y+i+(x+j)*512)[0] = group_number

            _tmp_ws = GroupDetectors(InputWorkspace=_tmp_ws, CopyGroupingFromWorkspace=_tmp_group, EnableLogging=False)
            DeleteWorkspace(_tmp_group, EnableLogging=False)

        progress.report('Adding logs')

        # Hack: ConvertToMD is needed so that a deep copy of the ExperimentInfo can happen
        # outWS.addExperimentInfo(_tmp_ws) # This doesn't work but should, when you delete `ws` `outWS` also loses it's ExperimentInfo
        _tmp_ws = Rebin(_tmp_ws, '0,1,2', EnableLogging=False)
        _tmp_ws = ConvertToMD(_tmp_ws, dEAnalysisMode='Elastic', EnableLogging=False, PreprocDetectorsWS='__PreprocessedDetectorsWS')

        preprocWS = mtd['__PreprocessedDetectorsWS']
        twotheta = preprocWS.column(2)
        azimuthal = preprocWS.column(3)

        outWS.copyExperimentInfos(_tmp_ws)
        DeleteWorkspace(_tmp_ws, EnableLogging=False)
        DeleteWorkspace('__PreprocessedDetectorsWS', EnableLogging=False)
        # end Hack

        add_time_series_property('s1', outWS.getExperimentInfo(0).run(), time_ns_array, s1_array)
        outWS.getExperimentInfo(0).run().getProperty('s1').units = 'deg'
        add_time_series_property('duration', outWS.getExperimentInfo(0).run(), time_ns_array, duration_array)
        outWS.getExperimentInfo(0).run().getProperty('duration').units = 'second'
        outWS.getExperimentInfo(0).run().addProperty('run_number', run_number_array, True)
        add_time_series_property('monitor_count', outWS.getExperimentInfo(0).run(), time_ns_array, monitor_count_array)
        outWS.getExperimentInfo(0).run().addProperty('twotheta', twotheta, True)
        outWS.getExperimentInfo(0).run().addProperty('azimuthal', azimuthal, True)

        setGoniometer_alg = self.createChildAlgorithm("SetGoniometer", enableLogging=False)
        setGoniometer_alg.setProperty("Workspace", outWS)
        setGoniometer_alg.setProperty("Axis0", 's1,0,1,0,1')
        setGoniometer_alg.setProperty("Average", False)
        setGoniometer_alg.execute()

        # Perform the normalization
        # TODO: Confirm the format of the normalization input (file or workspace or both?)
        # assuming we are getting a workspace
        ws_norm = self.getProperty("NormWorkspace")
        norm_method = self.getProperty("NormalizedBy").value

        if ws_norm and (norm_method != ''):
            print("Perform normalization")
            norm_replicated = ReplicateMD(ShapeWorkspace=outWS, DataWorkspace=ws_norm)
            outWS = DivideMD(LHSWorkspace=outWS, RHSWorkspace=norm_replicated)
            if norm_method.lower() == 'counts':
                scale = ws_norm.getSignalArray().mean()
                print('scale counts = {}'.format(int(scale)))
            elif norm_method.lower() == 'monitor':
                scale = np.array(outWS.getExperimentInfo(0).run().getProperty('monitor_count').value)
                scale /= ws_norm.getExperimentInfo(0).run().getProperty('monitor_count').value[0]
                print('scale monitor = {}'.format(scale))
            elif norm_method.lower() == 'time':
                scale = np.array(ws_norm.getExperimentInfo(0).run().getProperty('duration').value)
                scale /= ws_norm.getExperimentInfo(0).run().getProperty('duration').value[0]
            else:
                raise RuntimeError('Unknown normalization method: {}'.format(norm_method))
            # perform the normalization
            outWS.setSignalArray(outWS.getSignalArray() / scale)
            outWS.setErrorSquaredArray(outWS.getErrorSquaredArray() / scale**2)
        else:
            print("Skipping normalization")

        self.setProperty("OutputWorkspace", outWS)


def add_time_series_property(name, run, times, values):
    log = FloatTimeSeriesProperty(name)
    for t, v in zip(times, values):
        log.addValue(t, v)
    run[name] = log


AlgorithmFactory.subscribe(LoadWANDSCD)
