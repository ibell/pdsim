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
    cpdef resize(self, int Nvalues)
    cdef public arraym flows
    cdef public int Nflows
    cdef public __name__, __strf__
    
cdef class PyFlowFunctionWrapper(FlowFunction):
    cdef public dict kwargs
    cdef public object Function
    
    cpdef double call(self, FlowPath FP) except *
    
cdef class IsentropicNozzleWrapper(FlowFunction):
    cpdef double call(self, FlowPath FP) except *
    
cpdef double FrictionCorrectedIsentropicNozzle(double A, State State_up, State State_down, double delta, int Type, double t = *, double ro = *)

cdef class ValveModel(object):
    cdef public double A_port,A_valve,d_valve,rho_valve,d_port,m_eff,C_D,k_valve,x_stopper
    cdef public object key_up, key_down
    cdef public State State_up, State_down
    cdef public arraym xv
    cdef public double x_tr
    
    cpdef set_xv(self, arraym xv)
    cpdef arraym get_xv(self)
    cpdef double A(self)
    
    @cython.locals(exists_keys = cython.list)
    cpdef get_States(self, Core)
    
    cdef _pressure_dominant(self,arraym f, double x, double xdot, double rho, double V, double deltap)
    cdef _flux_dominant(self,arraym f, double x, double xdot, double rho, double V)
    cpdef double flow_velocity(self,State State_up, State State_down)
    cpdef arraym derivs(self, Core)
    cpdef dict __cdict__(self)
