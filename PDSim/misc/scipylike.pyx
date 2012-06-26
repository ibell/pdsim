import cython
import numpy as np
cimport numpy as np

def trapz(np.ndarray[np.float_t] y, np.ndarray[np.float_t] x):
    cdef int i,n
    cdef double sum=0
    
    n=y.size
    
    for i in range(n-1):
        sum += (y[i]+y[i+1])/2.0*(x[i+1]-x[i])
    
    return sum