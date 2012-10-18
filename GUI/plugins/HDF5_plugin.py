import sys
sys.path.append('..')

import pdsim_plugins
import pdsim_panels
import wx
import numpy as np
import types
import h5py

class StubClass():
    def __init__(self,d):
        self.__dict__.update(d)
        
class HDF5Writer(object):
    
    def _recursive_write(self, f, struct):
        
        for thing in dir(struct):
            #Skip everything starting with '_'
            if thing.startswith('_'):
                continue
            try:
                #Get the attribute
                value = getattr(struct, thing)
                
            except AttributeError:
                #If it can't get the attribute, just go to the next thing
                continue
            
            #If it is an integer, floating point value, or numpy array
            if isinstance(value,(int, float, np.ndarray)):
                #Save it as a value, go to next thing
                f.create_dataset(thing, data = value)
                continue
            elif isinstance(value, basestring):
                str_type = h5py.new_vlen(str)
                f.create_dataset(thing, dtype=str_type, data = value)
                continue
            
            import inspect
            #Skip methods, functions, built-in functions and routines
            if (inspect.ismethod(value)
                or inspect.isfunction(value)
                or inspect.isbuiltin(value)
                or inspect.isroutine(value)):
                    continue
            
            if type(value) is types.DictType:
                dict_group = f.create_group(thing)
                # Recurse into the entries in the dictionary by turning the 
                # dictionary into a class
                self._recursive_write(dict_group, StubClass(value))
            
            elif isinstance(value, list):
                dict_group = f.create_group(thing)
                #Convert to numpy array
                #List to a class
                cls = StubClass({str(i):v for i,v in enumerate(value)})
                #Write class recursively
                self._recursive_write(dict_group, cls)
            else:
                f.create_group(thing)
                #Recurse into the class
                self._recursive_write(f[thing], value)
        
    def write_to_file(self, struct, fName):
        """
        Write the structure to the file given by fName
        """
        
        f = h5py.File(fName,'w')
        self._recursive_write(f, struct)
        f.close()
    
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