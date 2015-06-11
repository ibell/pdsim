from CoolProp.State import State as StateClass
from CoolProp.State cimport State as StateClass

# For oil flooding
from PDSim.core.state_flooded import StateFlooded as StateClassFlood
from PDSim.core.state_flooded cimport StateFlooded as StateClassFlood

from libcpp cimport bool
from PDSim.core.core import PDSimCore
from PDSim.flow.flow import FlowPathCollection
from PDSim.flow.flow cimport FlowPathCollection
from PDSim.misc.datatypes cimport arraym
from PDSim.misc.datatypes import arraym

cpdef public bint __hasLiquid__ 

cdef class Tube(object):
    cdef readonly bool __hasLiquid__
    cdef public bytes key1,key2
    cdef public int fixed
    cdef public StateClass State1, State2
    
    # For oil flooding
    cdef public StateClassFlood StateFlood1, StateFlood2
    
    cdef public object TubeFcn
    cdef public double Q_add,alpha,L,ID,OD,mdot,Q
    cdef public bool exists
    cdef public int i1,i2

cdef class TubeCollection(list):
    cdef readonly bint __hasLiquid__
    cdef dict _Nodes
    cdef arraym harray, parray, Tarray, xLarray
    cpdef update_existence(self, int NCV)
    cpdef arraym get_h(self)
    cpdef arraym get_p(self)
    cpdef arraym get_T(self)
    cpdef arraym get_xL(self)
    cpdef dict get_Nodes(self)
    cpdef update(self)
    
cdef class CVScore(object):
    cdef readonly bint __hasLiquid__
    cdef list array_list
    cpdef update_size(self, int N)
    cdef build_all(self, int N)
    cdef free_all(self)
    cpdef copy(self)
    cpdef calculate_flows(self, FlowPathCollection Flows, arraym harray, arraym parray, arraym Tarray)
    cpdef calculate_flows_flood(self, FlowPathCollection Flows, arraym harray, arraym parray, arraym Tarray, arraym xLarray)

cdef class CVArrays(CVScore):
    cdef public arraym T,p,h,rho,V,dV,cp,cv,m,v,dpdT_constV,Q,xL,dudxL
    cdef public arraym drhodtheta, dTdtheta, dmdtheta, dxLdtheta, summerdm, summerdT, summerdxL, property_derivs
    cdef int state_vars, N
    cdef double omega
    cpdef just_volumes(self, list CVs, double theta)
    cpdef properties_and_volumes(self, list CVs, double theta, int state_vars, arraym x)
    cpdef calculate_derivs(self, double omega, bint has_liquid)
    
cdef class ControlVolume(object):
    cdef readonly bint __hasLiquid__
    cdef public long keyIndex
    cdef public bytes key, discharge_becomes
    cdef public object becomes
    cdef public object V_dV
    cdef public dict V_dV_kwargs
    cdef public object ForceFcn
    cdef public bint exists
    cdef public StateClass State
    
    # For oil flooding
    cdef public StateClassFlood StateFlood

cdef class ControlVolumeCollection(object):
    cdef readonly bint __hasLiquid__
    cdef readonly list keys, CVs, indices, exists_keys, exists_indices, exists_CV
    cdef readonly dict Nodes
    cdef readonly int N, Nexist
    cdef public CVArrays arrays
    cpdef add(self, ControlVolume CV)
    cpdef rebuild_exists(self)
    cpdef updateStates(self, str name1, arraym array1, str name2, arraym array2)
    cpdef updateStatesFlood(self, str name1, arraym array1, str name2, arraym array2, str name3, arraym array3)
    cpdef volumes(self, double theta, bint as_dict = *)
    cpdef at(self, int i)
