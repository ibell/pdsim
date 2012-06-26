
cdef double sqrt(double)
cdef double pow(double, double)

cdef class FlowPath:
    cdef public bytes key_up, key_down, key1, key2, Gas
    cdef public double p_up, p_down, T_up, h_up, A, mdot
    cdef public object MdotFcn, MdotFcn_kwargs

    @cython.locals(p1=cython.double, p2=cython.double)
    cpdef calculate(self)
    
cdef class FlowPathCollection(list):
    pass