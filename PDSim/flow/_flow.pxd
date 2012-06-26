from CoolProp.State import State
from CoolProp.State cimport State

cpdef tuple sum_flows(bytes key, object Flows)

cdef class _FlowPathCollection(list):
    pass
    
#Make a stripped down class with the necessary terms included
cdef class _FlowPath:
    cdef public double mdot, h_up, T_up, p_up, p_down, A
    cdef public bytes key_up, key_down, key1, key2, Gas
    cdef public dict MdotFcn_kwargs
    cdef public State State1,State2,State_up,State_down
    
    cpdef dict __cdict__(self,AddStates=*)
    cpdef _FlowPath get_deepcopy(self)
    cpdef _calculate(self)