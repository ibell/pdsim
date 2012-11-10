from CoolProp.State import State
from CoolProp.State cimport State

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.map cimport map
from cython.operator cimport dereference as deref, preincrement as inc
    
from PDSim.flow.flow_models import FlowFunctionWrapper
from PDSim.flow.flow_models cimport FlowFunctionWrapper

#cpdef sumterms_given_CV(bytes key, list flows)

cdef class FlowPathCollection(list):
    cpdef update_existence(self, Core)
    cpdef calculate(self, dict hdict)
    cpdef get_deepcopy(self)
    cpdef deepcopy(self)
    cpdef sumterms(self,Core)
        
#Make a stripped down class with the necessary terms included
cdef class FlowPath(object):
    cdef public bint exists, key1_exists, key2_exists, key_up_exists, key_down_exists
    cdef int ikey1, ikey2, ikey_up, ikey_down
    cdef public double mdot, h_up, T_up, p_up, p_down, A
    cdef public bytes key_up, key_down, key1, key2, Gas
    cdef public FlowFunctionWrapper MdotFcn
    cdef public State State1,State2,State_up,State_down
    
    cpdef dict __cdict__(self, AddStates = *)
    cpdef FlowPath get_deepcopy(self)
    cpdef calculate(self, dict hdict = *)