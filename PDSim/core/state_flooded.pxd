from libcpp.string cimport string
import cython
cimport cython

from libcpp.vector cimport vector

from CoolProp.CoolProp cimport State as State
from CoolProp.CoolProp cimport constants_header
#from CoolProp cimport AbstractState as cAbstractState
    
cdef extern from "Python.h":
    char* __FILE__

cdef extern from "Python.h":
    int __LINE__

cdef class StateFlooded(State):

    cdef readonly bytes Liq
    cdef readonly double xL_
    cdef bytes model
    #cdef int iLiquid    # Maybe ??
    
    cpdef speed_test(self, int N)
    cpdef update(self, dict params)
    cpdef update_ph(self, double p, double h)
    cpdef update_TrhoxL(self, double T, double rho, double xL)
    cpdef State copy(self)
    cpdef double Props(self, constants_header.parameters iOutput) except *
    cpdef double VoidFrac(self, string Ref, string Liq, double T, double P, double xL) except *
    cpdef long Phase(self) except *
    cpdef double get_T(self) except *
    cpdef double get_p(self) except *
    cpdef double get_Q_m(self) except *
    cpdef double get_h(self) except *
    cpdef double get_rho(self) except *
    cpdef double get_s(self) except *
    cpdef double get_u(self) except *
    cpdef double get_visc(self) except *
    cpdef double get_cond(self) except *
    cpdef double get_cp(self) except *
    cpdef double get_cp0(self) except *   #TODO: not sure if we need this one
    cpdef double get_cv(self) except *
    cpdef double get_dpdT(self) except *
    cpdef double get_dudxL(self) except*
    
    cpdef double k_mix(self, string, string, double, double, double) except *
    cpdef double mu_mix(self, string, string, double, double, double) except *
    cpdef double h_mix(self, string, string, double, double, double) except *
    cpdef double s_mix(self, string, string, double, double, double) except *
    cpdef double u_mix(self, string, string, double, double, double) except *
    cpdef double cp_mix(self, string, string, double, double, double) except *
    cpdef double cv_mix(self, string, string, double, double, double) except *
    cpdef double rho_mix(self, string, string, double, double, double) except *

    cpdef double k_liq(self, string, double) except *
    cpdef double mu_liq(self, string, double) except *
    cpdef double h_liq(self, string, double, double) except *
    cpdef double s_liq(self, string, double) except *
    cpdef double u_liq(self, string, double) except *
    cpdef double cp_liq(self, string, double) except *
    cpdef double rho_liq(self, string, double) except *