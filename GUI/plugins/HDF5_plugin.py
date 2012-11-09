import sys
sys.path.append('..')

import pdsim_plugins
import pdsim_panels
import wx
import numpy as np
import types
from PDSim.misc.hdf5 import HDF5Writer
    
class HDF5Plugin(pdsim_plugins.PDSimPlugin):
    
    short_description = "HDF5 file output for simulation data"
    
    def should_enable(self):
        """ Returns True if the plugin should be enabled """
        try:
            import h5py
            return True
        except ImportError:
            import warnings
            warnings.warn('HDF5Plugin was not loaded because the h5py package could not be imported')
            return False
    
    def activate(self, event = None):
        """ Activate the plugin """
        self._activated = not self._activated
        
    def apply(self):
        """
        Doesn't need to do anything at build time of simulation before it is run
        """
        pass