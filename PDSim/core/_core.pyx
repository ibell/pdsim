from __future__ import division

import cython
cimport cython

import numpy as np
cimport numpy as np

#Import the listm type used in PDSim
from PDSim.misc._listmath import listm
from PDSim.misc._listmath cimport listm

cpdef delimit_vector(np.ndarray[np.float_t, ndim = 1] x0, np.ndarray[np.float_t, ndim = 1] y0):
        """
        Break vectors into continuous chunks of real values as needed
        
        All the trailing NAN values have been removed already
        """
        cdef long i,L,R
        cdef list real_bounds,nan_bounds,x_list,y_list
        
        real_bounds = []
        nan_bounds = []
        
        if not np.any(np.isnan(y0)):
            return [x0],[y0]
        
        if np.isnan(y0[0]):
            nan_bounds += [0]
        else:
            real_bounds += [0]
            
        if np.isfinite(y0[-1]):
            real_bounds += [len(y0)-1]
        else:
            nan_bounds += [len(y0)-1]
            
        for i in range(len(y0)-1):
            if np.isnan(y0[i])==[False] and np.isnan(y0[i+1])==[True]:
                real_bounds+=[i]
                nan_bounds+=[i+1]
            elif np.isnan(y0[i])==[True] and np.isnan(y0[i+1])==[False]:
                nan_bounds += [i]
                real_bounds += [i+1]
        
        nan_bounds=sorted(nan_bounds)
        real_bounds=sorted(real_bounds)
        
        #Check that all the NAN chunks are all NAN in fact        
        for i in range(0,len(nan_bounds),2):
            L = nan_bounds[i]
            R = nan_bounds[i+1]+1
            if not np.all(np.isnan(y0[L:R])):
                raise ValueError('All the elements in NAN chunk are not NAN')
        
        x_list,y_list = [],[]
        #Calculate the 
        for i in range(0,len(real_bounds),2):
            L = real_bounds[i]
            R = real_bounds[i+1]+1
            x_list.append(x0[L:R])
            y_list.append(y0[L:R])
        
        return x_list, y_list
    
cdef class _PDSimCore:
    pass
    