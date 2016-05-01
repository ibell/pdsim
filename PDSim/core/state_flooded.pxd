from libcpp.string cimport string
import cython
cimport cython

from libcpp.vector cimport vector

from CoolProp.CoolProp cimport State as State
cimport CoolProp.constants_header as constants

cdef class StateFlooded(State):

    cdef readonly bytes Liq
    cdef readonly double xL_, alpha
    cdef bytes model
    #cdef int iLiquid    # Maybe ??
    
    cpdef speed_test(self, int N)
    cpdef update(self, dict params)
    cpdef update_ph(self, double p, double h)
    cpdef update_TrhoxL(self, double T, double rho, double xL)
    cpdef StateFlooded copy2(self)
    cpdef double Props(self, constants.parameters iOutput) except *
    
    cpdef long Phase(self) except *
    cpdef double get_T(self) except *
    cpdef double get_p(self) except *
    cpdef double get_Q(self) except *
   
    cpdef double get_hmix(self) except *
    cpdef double get_rhomix(self) except *
    cpdef double get_smix(self) except *
    cpdef double get_umix(self) except *
    cpdef double get_emix(self) except *
    cpdef double get_viscmix(self) except *
    cpdef double get_condmix(self) except *
    cpdef double get_Prmix(self) except *
    cpdef double get_cpmix(self) except *
    cpdef double get_cp0(self) except *   #TODO: not sure if we need this one
    cpdef double get_cvmix(self) except *
    cpdef double get_kstar(self) except *
    cpdef double get_VoidFrac(self) except *

    cpdef double get_hL(self) except *
    cpdef double get_rhoL(self) except *
    cpdef double get_sL(self) except *
    cpdef double get_uL(self) except *
    cpdef double get_viscL(self) except *
    cpdef double get_condL(self) except *
    cpdef double get_PrL(self) except *
    cpdef double get_cL(self) except *

    cpdef double get_dpdT(self) except *
    cpdef double get_dudxL(self) except *
    cpdef double get_cKe(self) except *
    cpdef double get_cve(self) except *
    
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
    cpdef double Pr_liq(self) except *
    cpdef double VoidFrac(self) except *
    
    cpdef double cK_e(self) except *
    cpdef double cv_e(self) except *
