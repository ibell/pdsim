from PDSim.misc.datatypes import arraym
from PDSim.misc.datatypes cimport arraym

cdef class _Recip(object):
    
    cdef public double crank_length, connecting_rod_length, A_piston, V_dead 
    cdef public double omega, piston_diameter
    
    #Function prototypes
    cpdef tuple V_dV(self, double theta)
    cpdef arraym heat_transfer_callback(self, double theta)
    cpdef dict __cdict__(self)