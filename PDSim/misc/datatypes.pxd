import cython
cimport cython
        
cdef class listm(list):
    """
    See http://docs.cython.org/src/userguide/special_methods.html
    """
    
cdef class arraym(object):

    cdef double* data
    cdef readonly int N
    
    cdef void set_data(self, double *data, int N)
    cpdef set_size(self, int N)
    cpdef arraym copy(self)