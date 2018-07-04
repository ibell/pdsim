from __future__ import print_function

import h5py
import numpy as np
import types
import six 

class StubClass():
    def __init__(self,d):
        self.__dict__.update(d)
        
class HDF5Writer(object):
    """
    This class contains the logic for writing a nested structure to a file.  The structure to be written could include dictionaries, lists, numpy arrays, classes, and combinations thereof.  The recursive writer function will walk through the tree, writing all the elements into the HDF5 file in the same general structure as was laid out in the object passed.
    
    Some small modifications are needed - for instance, lists are written in a slighltly different manner.  Otherwise, the structure should be pretty much exactly the same.
    """
    
    def _recursive_write(self, f, struct):
        
        for thing in dir(struct):
            # Skip everything starting with '_'
            if thing.startswith('_'):
                continue
            try:
                # Get the attribute
                value = getattr(struct, thing)
            except (AttributeError,ValueError) as E:
                print((thing, E))
                # If it can't get the attribute, just go to the next thing
                continue
            
            # If it is an integer, floating point value, or numpy array
            if isinstance(value, (int, float)):
                # Save it as a value, go to next thing
                f.create_dataset(thing, data = value)
                continue
            elif isinstance(value, np.ndarray):
                if not value.shape: # value.shape is an empty tuple
                    # It's a one-element numpy array, or empty
                    if not value or not value.shape or min(value.shape) == 0:
                        continue
                    else:
                        print(value, len(value), type(value))
                        f.create_dataset(thing, data = value)
                else:
                    # Save it with compression, go to next thing
                    f.create_dataset(thing, data = value, compression = 'gzip')
                continue
            elif isinstance(value, six.string_types):
                str_type = h5py.new_vlen(str)
                f.create_dataset(thing, dtype=str_type, data = value)
                continue
            
            import inspect
            # Skip methods, functions, built-in functions, routines, and modules
            if (inspect.ismethod(value)
                or inspect.isfunction(value)
                or inspect.isbuiltin(value)
                or inspect.isroutine(value)
                or inspect.ismodule(value)):
                    continue
            
            if isinstance(value, dict):
                dict_group = f.create_group(thing)
                # Recurse into the entries in the dictionary by turning the 
                # dictionary into a class
                self._recursive_write(dict_group, StubClass(value))
            
            elif isinstance(value, (list,tuple)):
                dict_group = f.create_group(thing)
                # Convert to numpy array
                # List/Tuple to a class
                cls = StubClass({str(i):v for i,v in enumerate(value)})
                # Write class recursively
                self._recursive_write(dict_group, cls)
            else:
                f.create_group(thing)
                # Recurse into the class
                self._recursive_write(f[thing], value)
        
    def write_to_file(self, struct, fName):
        """
        Write the structure to the file given by fName
        """
        
        f = h5py.File(fName,'w')
        self._recursive_write(f, struct)
        f.close()
        
    def prune(self, fName, keys):
        """
        Prune datasets and groups that are nested below the keys requested
        """
        
        f = h5py.File(fName,'r+') #open r/w
        for key in keys:
            if key in f:
                del f[key]
            else:
                print('Could not prune this key from HFD5 file:',key)
        f.close()
