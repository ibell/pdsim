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
    cdef readonly double xL_, alpha
    cdef bytes model
    #cdef int iLiquid    # Maybe ??
    
    cpdef speed_test(self, int N)
    cpdef update(self, dict params)
    cpdef update_ph(self, double p, double h)
    cpdef update_TrhoxL(self, double T, double rho, double xL)
    cpdef State copy(self)
    cpdef double Props(self, constants_header.parameters iOutput) except *
    cpdef double VoidFrac(self) except *
    cpdef long Phase(self) except *
    cpdef double get_T(self) except *
    cpdef double get_p(self) except *
    cpdef double get_Q(self) except *
    cpdef double get_h(self) except *
    cpdef double get_rho(self) except *
    cpdef double get_s(self) except *
    cpdef double get_u(self) except *
    cpdef double get_e(self) except*
    cpdef double get_visc(self) except *
    cpdef double get_cond(self) except *
    cpdef double get_cp(self) except *
    cpdef double get_cp0(self) except *   #TODO: not sure if we need this one
    cpdef double get_cv(self) except *
    cpdef double get_dpdT(self) except *
    cpdef double get_dudxL(self) except*
    cpdef double get_cKe(self) except*
    cpdef double get_cve(self) except*
    
    cpdef double s_mix(self) except *
    cpdef double u_mix(self) except *
    cpdef double h_mix(self) except *
    cpdef double e_mix(self) except *
    cpdef double rho_mix(self) except *
    cpdef double cp_mix(self) except *
    cpdef double cv_mix(self) except *
    cpdef double mu_mix(self) except *
    cpdef double k_mix(self) except * 
    cpdef double Pr_mix(self) except * 
    cpdef double kstar_mix(self) except * 
    cpdef double dudxL_mix(self) except *
    cpdef double dpdT_const_V(self) except *
    cpdef double T_crit(self) except *
    cpdef double p_crit(self) except *

    cpdef double s_liq(self) except *
    cpdef double u_liq(self) except *
    cpdef double h_liq(self) except *
    cpdef double rho_liq(self) except *
    cpdef double cp_liq(self) except *
    cpdef double mu_liq(self) except *
    cpdef double k_liq(self) except *
    
    cpdef double cK_e(self) except *
    cpdef double cv_e(self) except *
