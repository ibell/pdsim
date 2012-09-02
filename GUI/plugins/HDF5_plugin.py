import sys
sys.path.append('..')

import pdsim_plugins
import pdsim_panels
import wx
import numpy as np

class HDF5Writer(object):
    
    def _recursive_write(self,f,struct):
        
        for thing in dir(struct):
            #Skip everything starting with '_'
            if thing.startswith('_'):
                continue
            
            try:
                #Get the attribute
                value = getattr(struct, thing)
            except AttributeError:
                print "couldn't get",thing
                #If it can't get the attribute, just go to the next thing
                continue
            
            #If it is an integer, floating point value, or numpy array
            if isinstance(value,(int, float, np.ndarray)):
                #Save it as a value, go to next thing
                f.create_dataset(thing, data = value)
                continue
            elif isinstance(value,list):
                #Convert to numpy array
                print np.array(value)
                f.create_dataset(thing, data = np.array(value))
            
            import inspect
            #Skip methods, functions, built-in functions and routines
            if (inspect.ismethod(value)
                or inspect.isfunction(value)
                or inspect.isbuiltin(value)
                or inspect.isroutine(value)):
                    continue
            
            import types
            print thing, type(value),type(value) is types.InstanceType
            
            if type(value) is types.DictType:
                f.create_group(thing)
                for k,v in value.iteritems():
                    f[thing].create_dataset(k, data = v)
            else:
                f.create_group(thing)
                #Recurse into the class
                self._recursive_write(f[thing], value)
        
    def write_to_file(self, struct, fName):
        """
        Write the structure to the file given by fName
        """
        import h5py
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