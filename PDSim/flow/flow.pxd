from CoolProp.State import State
from CoolProp.State cimport State

# For oil flooding
from PDSim.core.state_flooded import StateFlooded
from PDSim.core.state_flooded cimport StateFlooded

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.map cimport map
from cython.operator cimport dereference as deref, preincrement as inc
    
from PDSim.flow.flow_models import FlowFunction
from PDSim.flow.flow_models cimport FlowFunction

from PDSim.misc.datatypes import arraym
from PDSim.misc.datatypes cimport arraym

cdef class FlowPathCollection(list):
    cdef readonly bint __hasLiquid__
    cdef int N
    cdef int Nexists
    
    # The rotational speed [rad/s]
    cdef double omega
    
    cpdef update_existence(self, Core)
    cpdef calculate(self, arraym harray, arraym parray, arraym Tarray)
    cpdef calculateFlood(self, arraym harray, arraym parray, arraym Tarray, arraym xLarray)
    cpdef get_deepcopy(self)
    cpdef sumterms(self, arraym summerdT, arraym summerdm)
    cpdef sumtermsFlood(self, arraym summerdT, arraym summerdm, arraym summerdxL)
    cpdef list flow_paths
        
#Make a stripped down class with the necessary terms included
cdef class FlowPath(object):
    cdef public bytes key_up
    """The string key corresponding to the upstream node"""
    
    cdef public bytes key_down
    """ The string key corresponding to the downstream node """
    
    cdef public bytes key1
    """ The string key corresponding to the first node """
    
    cdef public bytes key2
    """ The string key corresponding to the second node """
    
    cdef public bint exists, key1_exists, key2_exists, key_up_exists, key_down_exists
    cdef public long key1Index, key2Index, key_up_Index, key_down_Index
    cdef public int ikey1, ikey2, ikey_up, ikey_down
    
    cdef public double mdot
    """ The mass flow rate [kg/s]"""
    
    cdef public double xL
    
    cdef public double h_up
    """ The upstream enthalpy [kJ/kg] """ 
    
    cdef public double h_down
    """ The downstream enthalpy [kJ/kg] """
    
    cdef public double T_up
    """ The upstream temperature [K] """
    
    cdef public double p_up
    """ The upstream pressure [kPa] """
    
    cdef public double p_down
    """ The downstream pressure [kPa] """
    
    cdef public double A
    """ The flow area [m^2] """
    
    cdef public double edot
    """ The rate of irreversibility generation in this flow path [kW]"""
    
    cdef public FlowFunction MdotFcn
    """ The function that will return the mass flow rate """
    
    cdef public bytes MdotFcn_str
    
    cdef public State State1
    """ The first state """
    
    cdef public State State2
    """ The second state """
    
    cdef public State State_up
    """ The upstream state """
    
    cdef public State State_down
    """ The downstream state """
    
    cpdef dict __cdict__(self, AddStates = *)
    cpdef FlowPath get_deepcopy(self)
    cpdef calculate(self, arraym harray, arraym parray, arraym Tarray)
    cpdef calculateFlood(self, arraym harray, arraym parray, arraym Tarray, arraym xLarray)
