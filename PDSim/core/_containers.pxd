from CoolProp.State import State as StateClass
from CoolProp.State cimport State as StateClass

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
    cdef list array_list
    
    #Storage arrays
    cdef public arraym T,p,h,rho,V,dV,cp,cv,m,v,dpdT_constV,Q,xL,dudxL
    
    #Property derivative arrays
    cdef public arraym drhodtheta, dTdtheta, dmdtheta, dxLdtheta, summerdm, summerdT, summerdxL, property_derivs
    
    # Other variables
    cdef int state_vars,N
    cdef double omega
    
    cpdef update_size(self, int N)
    cdef build_all(self, int N)
    cdef free_all(self)
    
    cpdef properties_and_volumes(self, CVs, double theta, int state_vars, arraym x)
    #
    cpdef calculate_flows(self, Flows, harray, Core)
    
    cpdef calculate_derivs(self, double omega, bint has_liquid)
    
    cpdef copy(self)
    