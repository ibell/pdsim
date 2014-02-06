import cython
cimport cython

import numpy as np
cimport numpy as np

from libcpp.vector cimport vector

from CoolProp.State cimport State

from PDSim.flow.flow import FlowPath
from PDSim.flow.flow cimport FlowPath

from PDSim.misc.datatypes import arraym, empty_arraym
from PDSim.misc.datatypes cimport arraym, empty_arraym
    
from libc.math cimport exp, log, M_PI as pi, M_E as e, sqrt

cpdef double IsentropicNozzle(double A, State State_up, State State_down, int other_output=*)

cdef class FlowFunction(object):
    cpdef double call(self, FlowPath FP) except *
    
cdef class PyFlowFunctionWrapper(FlowFunction):
    cdef dict kwargs
    cdef public object Function
    
    cpdef double call(self, FlowPath FP) except *
    
cdef class IsentropicNozzleWrapper(FlowFunction):
    cpdef double call(self, FlowPath FP)
    
cpdef double FrictionCorrectedIsentropicNozzle(double A, State State_up, State State_down, double delta, int Type, double t = *, double ro = *)

cdef class ValveModel(object):
    cdef public double E,A_port,A_valve,d_valve,l_valve,a_valve,h_valve,rho_valve,d_port,m_eff,C_D,k_valve,x_stopper
    cdef public bytes key_up, key_down
    cdef public State State_up, State_down
    cdef public arraym xv
    cdef public double x_tr
    
    cpdef set_xv(self, arraym xv)
    cpdef double A(self)
    
    @cython.locals(exists_keys = cython.list, key = cython.bytes)
    cpdef get_States(self, Core)
    
    cdef _pressure_dominant(self,arraym f, double x, double xdot, double rho, double V, double deltap)
    cdef _flux_dominant(self,arraym f, double x, double xdot, double rho, double V)
    cpdef double flow_velocity(self,State State_up, State State_down)
    cpdef arraym derivs(self, Core)
    cpdef dict __cdict__(self)
