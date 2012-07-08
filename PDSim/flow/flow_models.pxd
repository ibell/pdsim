import cython
cimport cython

import numpy as np
cimport numpy as np

from CoolProp.State cimport State

@cython.locals(
cp=cython.double,
cv=cython.double,
A=cython.double, 
T_up=cython.double,
p_up=cython.double, 
p_down=cython.double, 
k=cython.double,
R=cython.double,
rho_up=cython.double,
pr=cython.double,
pr_crit=cython.double,
v=cython.double,
c=cython.double,
e=cython.double,
f=cython.double,
mdot=cython.double,
T_down=cython.double,
rho_down=cython.double,
otherparameters=cython.dict,)
cpdef IsentropicNozzle(double A, State State_up, State State_down, bint full_output=*)

from libc.math cimport exp, log, M_PI as pi, M_E as e, sqrt

@cython.locals(Re = cython.double, v = cython.double)
cpdef FrictionCorrectedIsentropicNozzle(double A, State State_up, State State_down, double delta, bytes Type, double t = *, double ro = *, bint full_output = *)

cdef class ValveModel(object):
    cdef public double E,A_port,A_valve,d_valve,l_valve,a_valve,h_valve,rho_valve,d_port,m_eff,C_D,k_valve,x_stopper
    cdef public bytes key_up, key_down
    cdef public State State_up, State_down
    cdef public list xv
    
    cpdef set_xv(self, list xv)
    cpdef double A(self)
    
    @cython.locals(exists_keys = cython.list, key = cython.bytes)
    cpdef get_States(self, Core)
    
    cpdef tuple _pressure_dominant(self,double x, double xdot, double rho, double V, double deltap)
    cpdef tuple _flux_dominant(self,double x, double xdot, double rho, double V)
    cpdef double flow_velocity(self,State State_up, State State_down)
    
    @cython.locals(V=cython.double,x=cython.double,xdot=cython.double,x_tr=cython.double)
    cpdef list derivs(self, Core)
    cpdef dict __cdict__(self)


