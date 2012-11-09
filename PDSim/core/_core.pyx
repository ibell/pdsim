from __future__ import division

import cython
cimport cython

import numpy as np
cimport numpy as np

#Import the listm type used in PDSim
from PDSim.misc._listmath import listm
from PDSim.misc._listmath cimport listm

cpdef delimit_vector(np.ndarray[np.float_t, ndim = 1] x, np.ndarray[np.float_t, ndim = 1] y):
        """
        Break vectors into continuous chunks of real values as needed by splitting based on the y numpy array
        
        Parameters
        ----------
        x : numpy array
          May not contain NAN values  
        y : numpy array
        """
        
        #Return the arrays back if there are no NAN values in y
        if not np.any(np.isnan(y)):
            return [x],[y]
        
        #A list of all the transitions between NAN and non-NAN values in y vector
        #Values are the indices at the left
        indices = list(np.where(np.diff(np.isnan(y)))[0])
        
        segments = []
        
        if not np.isnan(y[0]):
            segments.append((0,indices[0]))
            keep = False
        else:
            keep = True
            
        for i in range(len(indices)-1):
            if keep:
                segments.append((indices[i]+1,indices[i+1]))
            
            #Flip the keep flag
            keep = not keep
        
        if indices[-1] < len(y)-1:
            if keep:
                segments.append((indices[-1]+1,len(y)-1))
            
        return [x[L:R+1] for L,R in segments], [y[L:R+1] for L,R in segments]
    
cdef class _PDSimCore:
    pass
    