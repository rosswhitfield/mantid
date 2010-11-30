#!/usr/bin/env python
import os

from distutils.core import setup

# Compile resource file
try:
    os.system("pyrcc4 -o reduction_gui/settings/qrc_resources.py reduction_gui/settings/resources.qrc")
    os.system("pyuic4 -o ui/ui_reduction_main.py ui/reduction_main.ui")
    os.system("pyuic4 -o ui/ui_hfir_background.py ui/hfir_background.ui")
    os.system("pyuic4 -o ui/ui_hfir_beam_finder.py ui/hfir_beam_finder.ui")
    os.system("pyuic4 -o ui/ui_hfir_data_summary.py ui/hfir_data_summary.ui")
    os.system("pyuic4 -o ui/ui_hfir_output.py ui/hfir_output.ui")
    os.system("pyuic4 -o ui/ui_hfir_transmission.py ui/hfir_transmission.ui")
    os.system("pyuic4 -o ui/ui_trans_direct_beam.py ui/trans_direct_beam.ui")
    os.system("pyuic4 -o ui/ui_trans_spreader.py ui/trans_spreader.ui")
    os.system("pyuic4 -o ui/ui_instrument_dialog.py ui/instrument_dialog.ui")
except:
    print "Could not compile resource file"
    
setup(name='SANSReduction',
      version='1.0',
      description='SANS Reduction for Mantid',
      author='Mathieu Doucet',
      author_email='doucetm@ornl.gov',
      url='http://www.mantidproject.org',
      packages=['reduction_gui', 'reduction_gui.widgets', 'reduction_gui.instruments', 
                'reduction_gui.reduction', 'reduction_gui.settings'],
      package_dir={'reduction_gui' : '',
                   'reduction_gui.widgets' : 'reduction_gui/widgets',
                   'reduction_gui.instruments' : 'reduction_gui/instruments',
                   'reduction_gui.settings' : 'reduction_gui/settings',
                   'reduction_gui.reduction' : 'reduction_gui/reduction'}
     )
