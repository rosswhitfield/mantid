from __future__ import (absolute_import, division, print_function)
from mantid.api import PythonAlgorithm, AlgorithmFactory, PropertyMode, WorkspaceProperty, FileProperty, FileAction
from mantid.kernel import Direction, IntArrayProperty, FloatTimeSeriesProperty
import numpy as np
import datetime
import os


class HB2AReduce(PythonAlgorithm):
    _gaps = np.array([0.   ,    2.641,    5.287,    8.042,   10.775,   13.488,
                      16.129,   18.814,   21.551,   24.236,   26.988,   29.616,
                      32.312,   34.956,   37.749,   40.4  ,   43.111,   45.839,
                      48.542,   51.207,   53.938,   56.62 ,   59.286,   61.994,
                      64.651,   67.352,   70.11 ,   72.765,   75.492,   78.204,
                      80.917,   83.563,   86.279,   88.929,   91.657,   94.326,
                      97.074,   99.784,  102.494,  105.174,  107.813,  110.551,
                      113.25 ,  115.915])

    def category(self):
        return 'DataHandling\\Nexus'

    def seeAlso(self):
        return [ "" ]

    def name(self):
        return 'HB2AReduce'

    def summary(self):
        return 'Performs data reduction for HB-2A POWDER at HFIR'

    def PyInit(self):
        self.declareProperty(FileProperty(name="Filename", defaultValue="",action=FileAction.Load,extensions=[".dat"]), "File to load")
        self.declareProperty(FileProperty(name="Vanadium",defaultValue="",action=FileAction.OptionalLoad,
                                          extensions=[".dat", ".txt"]),
                             doc="Vanadium file, can be either the vanadium scan file or the reduced vcorr file. "
                             "If not provided the vcorr file adjacent to the data file will be used")
        self.declareProperty(FileProperty(name="GapsFile",defaultValue="",action=FileAction.OptionalLoad,
                                          extensions=[".txt"]),
                             doc="If provided will use this gaps file otherwise uses default values")
        self.declareProperty(IntArrayProperty("ExcludeDetectors", []),
                             doc="Detectors to exclude. If not provided the HB2A_?__exclude_detectors.txt adjacent "
                             "to the data file will be used if it exist")
        self.declareProperty(WorkspaceProperty("IntermediateWorkspace", "",
                                               optional=PropertyMode.Optional,
                                               direction=Direction.Output),
                             "This workspace will include each anode as a separate spectrum, useful for debugging issues")
        self.declareProperty('Scale', 1.0, "The output will be scaled by this")
        self.declareProperty('BinData', False, "Data will be binned using BinSize. If False then all data will be unbinned")
        self.declareProperty('BinSize', 0.0505, "Bin size of the output workspace")
        self.declareProperty('OldStyle', False, "Use old style of binning, the same way that Graffiti does it")
        self.declareProperty(WorkspaceProperty("OutputWorkspace", "",
                                               optional=PropertyMode.Mandatory,
                                               direction=Direction.Output),
                             "Output Workspace")

    def validateInputs(self):
        issues = dict()

        return issues

    def PyExec(self):
        filename = self.getProperty("Filename").value
        scale = self.getProperty("Scale").value

        regexp = '#(.*)=(.*)'
        metadata = dict(np.char.strip(np.fromregex(filename, regexp, dtype='str')))

        indir, data_filename = os.path.split(filename)
        _, exp, scan = data_filename.replace(".dat", "").split('_')

        if self.getProperty("GapsFile").value:
            self._gaps = np.cumsum(np.genfromtxt(self.getProperty("GapsFile").value, usecols=(0)))

        detector_mask = np.full(44, True)
        if len(self.getProperty("ExcludeDetectors").value) == 0:
            exclude_filename = 'HB2A_{}__exclude_detectors.txt'.format(exp)
            exclude_detectors = np.loadtxt(os.path.join(indir, exclude_filename), ndmin=1, dtype=int)
        else:
            exclude_detectors = np.array(self.getProperty("ExcludeDetectors").value)
        detector_mask[exclude_detectors] = False

        data = np.genfromtxt(filename, skip_header=28, names=True)
        counts = np.array([data['anode{}'.format(n)] for n in range(1,45)])[detector_mask]

        twotheta = data['2theta']

        monitor = data['monitor']

        vanadium_count, vanadium_monitor, vcorr = self.get_vanadium(detector_mask,
                                                                    data['m1'][0],
                                                                    data['colltrans'][0],
                                                                    exp,
                                                                    indir)

        x = twotheta+self._gaps[:, np.newaxis][detector_mask]

        if self.getPropertyValue("IntermediateWorkspace"):
            # Separate spectrum per anode
            y, e = self.process(counts, scale, monitor, vanadium_count, vanadium_monitor, vcorr)
            createWS_alg = self.createChildAlgorithm("CreateWorkspace")
            createWS_alg.setProperty("DataX",x)
            createWS_alg.setProperty("DataY",y)
            createWS_alg.setProperty("DataE",e)
            createWS_alg.setProperty("NSpec",len(x))
            createWS_alg.setProperty("UnitX","Degrees")
            createWS_alg.setProperty("YUnitLabel","Counts")
            createWS_alg.setProperty("WorkspaceTitle",str(metadata['scan_title']))
            createWS_alg.execute()
            anodeOutWS = createWS_alg.getProperty("OutputWorkspace").value
            self.setProperty("IntermediateWorkspace", anodeOutWS)
            self.add_metadata(anodeOutWS, metadata, data)

        if self.getProperty("BinData").value:
            # Data binned with bin
            x, y, e = self.process_binned(counts, x.ravel(), scale, monitor, vanadium_count, vanadium_monitor, vcorr)
        else:
            # Sort array, all in one spectrum
            index_array = np.argsort(x.ravel())
            x, y, e = self.process(counts, scale, monitor, vanadium_count, vanadium_monitor, vcorr)
            x=x.ravel()[index_array]
            y=y.ravel()[index_array]
            e=e.ravel()[index_array]

        createWS2_alg = self.createChildAlgorithm("CreateWorkspace")
        createWS2_alg.setProperty("DataX",x)
        createWS2_alg.setProperty("DataY",y)
        createWS2_alg.setProperty("DataE",e)
        createWS2_alg.setProperty("UnitX","Degrees")
        createWS2_alg.setProperty("YUnitLabel","Counts")
        createWS2_alg.setProperty("WorkspaceTitle",str(metadata['scan_title']))
        createWS2_alg.execute()
        outWS = createWS2_alg.getProperty("OutputWorkspace").value

        self.setProperty("OutputWorkspace", outWS)

        self.add_metadata(outWS, metadata, data)

    def get_vanadium(self, detector_mask, m1, colltrans, exp, indir):
        """
        This function return either (vanadium_count, vanadium_monitor, None) or
        (None, None, vcorr) depending what type of file is provided by getProperty("Vanadium")
        """
        vanadium_filename = self.getProperty("Vanadium").value
        if vanadium_filename:
            if vanadium_filename.split('.')[-1] == 'dat':
                vanadium = np.genfromtxt(vanadium_filename)
                vanadium_count = vanadium[:, 5:49].sum(axis=0)[detector_mask]
                vanadium_monitor = vanadium[:, 3].sum()
                return vanadium_count, vanadium_monitor, None
            else:
                vcorr_filename = vanadium_filename
        else:
            # m1 = 0 -> Ge 115, 1.54A
            # m1 = 9.45 -> Ge 113, 2.41A
            # colltrans = 0 -> IN
            # colltrans = +/-80 -> OUT
            vcorr_filename = 'HB2A_{}__Ge_{}_{}_vcorr.txt'.format(exp,
                                                                  115 if np.isclose(m1, 0, atol=0.1) else 113,
                                                                  "IN" if np.isclose(colltrans, 0, atol=0.1) else "OUT")

        return None, None, np.genfromtxt(os.path.join(indir, vcorr_filename))[detector_mask]

    def process(self, counts, scale, monitor, vanadium_count=None, vanadium_monitor=None, vcorr=None):
        if vcorr is not None:
            y = counts/vcorr[:, np.newaxis]*scale/monitor
            e = np.sqrt(counts)/vcorr[:, np.newaxis]*scale/monitor
        else:
            y = counts/vanadium_count[:, np.newaxis]*vanadium_monitor/monitor
            e = np.sqrt(1/counts + 1/vanadium_count[:, np.newaxis] + 1/vanadium_monitor + 1/monitor)*y
        return y, e

    def process_binned(self, counts, x, scale, monitor, vanadium_count=None, vanadium_monitor=None, vcorr=None):
        binSize = self.getProperty("BinSize").value
        old = bool(self.getProperty("OldStyle").value)
        if old:
            y, _ = self.process(counts, 1.0, monitor, vanadium_count, vanadium_monitor, vcorr)
            y = y.flatten()
            x_out = []
            y_out = []
            e_out = []

            while len(x)>0:
                x_store = x.min() # Get smallest x value
                mask = np.logical_and(x >= x_store, x <= x_store+binSize) # Find every x value within binSize
                x_out.append(np.mean(x[mask])) # Append average y
                y_out.append(np.mean(y[mask])) # Append average x
                e_out.append(np.sqrt(np.sum(y[mask]))/sum(mask)) # Append average error
                x = x[np.logical_not(mask)] # Remove used values
                y = y[np.logical_not(mask)] # Remove used values
            x=np.array(x_out)
            y=np.array(y_out)
            e=np.array(e_out)
        else:
            bins = np.arange(x.min(), x.max()+binSize, binSize) # calculate bin boundaries
            inds = np.digitize(x, bins) # get bin indices

            # because np.broadcast_to is not in numpy 1.7.1 we use stride_tricks
            if vcorr is not None:
                vcorr=np.lib.stride_tricks.as_strided(vcorr, shape=counts.shape, strides=(vcorr.strides[0],0))
                vcorr_binned = np.bincount(inds, weights=vcorr.ravel(), minlength=len(bins))
            else:
                vanadium_count=np.lib.stride_tricks.as_strided(vanadium_count, shape=counts.shape, strides=(vanadium_count.strides[0],0))
                vanadium_binned = np.bincount(inds, weights=vanadium_count.ravel(), minlength=len(bins))
                vanadium_monitor_binned = np.bincount(inds, minlength=len(bins))*vanadium_monitor

            monitor=np.lib.stride_tricks.as_strided(monitor, shape=counts.shape, strides=(monitor.strides[0],0))

            counts_binned = np.bincount(inds, weights=counts.ravel(), minlength=len(bins))
            monitor_binned = np.bincount(inds, weights=monitor.ravel(), minlength=len(bins))
            number_binned = np.bincount(inds, minlength=len(bins))

            if vcorr is not None:
                y = np.nan_to_num(counts_binned/vcorr_binned*number_binned/monitor_binned)[1:]
                e = np.nan_to_num(np.sqrt(1/counts_binned)[1:])*y
            else:
                y = np.nan_to_num(counts_binned/vanadium_binned*vanadium_monitor_binned/monitor_binned)[1:]
                e = np.nan_to_num(np.sqrt(1/counts_binned + 1/vanadium_binned + 1/vanadium_monitor + 1/monitor_binned)[1:])*y
            x = bins
        return x, y*scale, e*scale

    def add_metadata(self, ws, metadata, data):
        run = ws.getRun()

        # Just copy all metadata in file
        for key in metadata.keys():
            if key != 'time': # time is used below as a TimeSeriesLog
                run.addProperty(key, str(metadata[key]), True)

        # Add correct start and end time
        start_time = datetime.datetime.strptime(metadata['time'] + ' ' + metadata['date'], '%I:%M:%S %p %m/%d/%Y')
        run.addProperty('start_time', start_time.isoformat(), True)

        time = data['time']
        time_array = [start_time + datetime.timedelta(seconds=s) for s in np.cumsum(time)]
        run.addProperty('end_time', time_array[-1].isoformat(), True)
        run.addProperty('duration', (time_array[-1]-time_array[0]).total_seconds(), True)

        # Add time series logs
        for name in data.dtype.names:
            if 'anode' not in name:
                log = FloatTimeSeriesProperty(name)
                for t, v in zip(time_array, data[name]):
                    log.addValue(t.isoformat(), v)
                run[name]=log


AlgorithmFactory.subscribe(HB2AReduce)
