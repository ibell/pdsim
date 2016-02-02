import cython
cimport cython

from cpython cimport bool
        
cdef class AnnotatedValue(object):
    cdef public object value
    cdef public str annotation, units, key
    
    
cdef class Collector(object):
    cdef public list vec
    """ The list that contains the values """
    
    cpdef v(self, int ndmin = *)
    
cdef class listm(list):
    """
    See http://docs.cython.org/src/userguide/special_methods.html
    """

cdef class arraym(object):

    cdef double* data
    cdef readonly int N
    """ The number of entries in this arraym instance """
    
    cdef void set_data(self, double *data, int N)
    cpdef set_size(self, int N)
    cpdef dealloc(self)
    cpdef arraym copy(self)
    cdef arraym slice(self, int i, int j)
    cpdef extend(self, arraym array2)
    cpdef double get_index(self, int i) except *
    cpdef double set_index(self, int i, double val) except *
    cpdef fill(self, double fillval)
    cpdef bool all_finite(self)
    
cpdef arraym empty_arraym(int N)
    
cdef inline check_dims(arraym x, arraym y):
    if x.N != y.N:
        raise ValueError('Cannot apply unary operator to arraym instances with lengths of '+str(x.N)+' and '+str(y.N))