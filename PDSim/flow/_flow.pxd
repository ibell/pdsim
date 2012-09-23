from CoolProp.State import State
from CoolProp.State cimport State
    
from PDSim.flow.flow_models import FlowFunctionWrapper
from PDSim.flow.flow_models cimport FlowFunctionWrapper

cpdef sumterms_given_CV(bytes key, list flows)

cdef class FlowPathCollection(list):
    cpdef calculate(self, Core, dict hdict)
    cpdef get_deepcopy(self)
    cpdef deepcopy(self)
    cpdef sumterms_helper(self, list exists_keys, double omega)
    cpdef sumterms(self,Core)
        
#Make a stripped down class with the necessary terms included
cdef class FlowPath(object):
    cdef public double mdot, h_up, T_up, p_up, p_down, A
    cdef public bytes key_up, key_down, key1, key2, Gas
    cdef public FlowFunctionWrapper MdotFcn
    cdef public State State1,State2,State_up,State_down
    
    cpdef dict __cdict__(self, AddStates = *)
    cpdef FlowPath get_deepcopy(self)
    cpdef calculate(self, dict hdict = *)