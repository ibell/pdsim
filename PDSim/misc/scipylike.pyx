import cython
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector

def trapz(np.ndarray[np.float_t] y, np.ndarray[np.float_t] x):
    """
    Do a trapezoidal integration of the array y with respect to x 
    
    Parameters
    ----------
    x : numpy array (1-D, float type)
    y : numpy array (1-D, float type)
    """
    cdef int i,n
    cdef double sum=0
    
    n=y.size
    
    for i in range(n-1):
        sum += (y[i]+y[i+1])/2.0*(x[i+1]-x[i])
    
    return sum
    
# A header defining the Spline class
cimport cSpline 
    
ctypedef fused number_t:
    double
    vector[double]
    
cdef class Spline:
    """ A python wrapper around the Spline interpolator class from Devin Lane: http://shiftedbits.org/2011/01/30/cubic-spline-interpolation/ """
    cdef cSpline.Spline[double,double] *thisptr     # hold a C++ instance which we're wrapping
    
    def __cinit__(self, vector[double] x, vector[double] y):
        self.thisptr = new cSpline.Spline[double, double](x, y)
        
    def __dealloc__(self):
        del self.thisptr
        
    def interpolate(self, number_t x):
        # Only one of these branches will be compiled for each specialization!
        if number_t is double:
            return self.thisptr.interpolate(x)
        else:
            # Select the vector<double> specialization of the interpolation function
            return np.array(self.thisptr.interpolate_vec( x ))
        
cpdef Spline splrep(vector[double] x, vector[double] y, k = 3, s = 0):
    """
    A C++/python equivalent of the scipy function scipy.interpolate.splrep
    """
    return Spline(x, y)
    

    
cpdef splev(number_t x, Spline spl):
    """
    A C++/python equivalent of the scipy function scipy.interpolate.splev
    """
    
    # Only one of these branches will be compiled for each specialization!
    if number_t is double:
        return spl.interpolate(x)
    else:
        return np.array(spl.interpolate_vec(x))
