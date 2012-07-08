from PDSim.misc._listmath import listm
from PDSim.misc._listmath cimport listm

cdef class _Recip(object):
    
    cdef public double crank_length, connecting_rod_length, A_piston, V_dead 
    cdef public double omega, piston_diameter
    
    #Function prototypes
    cpdef tuple V_dV(self, double theta)
    cpdef listm heat_transfer_callback(self, double theta)
    cpdef dict __cdict__(self)