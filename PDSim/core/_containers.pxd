from CoolProp.State cimport State as StateClass
from PDSim.misc._listmath cimport listm
from libcpp cimport bool

from PDSim.misc.datatypes cimport arraym
from PDSim.misc.datatypes import arraym
    
cdef class _ControlVolume:
    cdef public StateClass State
    
cdef class TubeCollection(list):
    cdef dict _Nodes
    cdef arraym harray
    
    cpdef update_existence(self, int NCV)
    cpdef arraym get_h(self)
    cpdef dict get_Nodes(self)
    cpdef update(self)
    
cdef class _Tube(object):
    cdef public bytes key1,key2
    cdef public int fixed
    cdef public double Q_add,alpha,L,ID,OD,mdot
    cdef public bool exists
    
cdef class CVArrays(object):
    #Storage arrays
    cdef public arraym T,p,h,rho,V,dV,cp,cv,m,v,dpdT_constV,Q,xL,dudxL
    
    # Other variables
    cdef int state_vars,N
    cdef double omega
    
    #Property derivative arrays
    cdef public arraym drhodtheta, dTdtheta, dmdtheta, dxLdtheta, summerdm, summerdT, summerdxL, property_derivs
    
    #
    cpdef calculate_flows(self, Flows, harray, Core)
    
    cpdef calculate_derivs(self, double omega, bint has_liquid)
    