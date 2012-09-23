#cython: embedsignature=True

import cython
cimport cython

import numpy as np
cimport numpy as np

from PDSim.misc._listmath import listm
from PDSim.misc._listmath cimport listm

cpdef list getcol(np.ndarray[np.float_t, ndim=2] mat, int colIndex, list indices):
    cdef int i,j
    cdef list values = []
    for i in indices: 
        values.append(mat[i,colIndex])
    #print 'getcol,i=',i,'indices=',indices,'vals=',values
    return values

cpdef setcol(np.ndarray[np.float_t, ndim=2] mat, int colIndex, list indices, listm x):
    cdef int i,j
    cdef double val
    #print 'setcol,i=',colIndex,'indices=',indices,'vals=',x
    if not len(indices) == len(x):
        print 'indices:',indices, 'x:',x
        raise ValueError('Length of indices [{0:d}] and x [{1:d}] are not the same'.format(len(indices),len(x)))
    for j in range(len(x)):
        i = indices[j]
        val = x[j]
        mat[i,colIndex]=val
    return None
