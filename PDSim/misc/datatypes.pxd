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
    cdef arraym slice(self, int i, int j)
    cpdef extend(self, arraym array2)
    
cdef inline check_dims(arraym x, arraym y):
    if x.N != y.N:
        raise ValueError('Cannot apply unary operator to arraym instances with lengths of '+str(x.N)+' and '+str(y.N))