from libcpp.string cimport string
import cython
cimport cython

from libcpp.vector cimport vector

include "AbstractState.pxd"
    
cdef extern from "Python.h":
    char* __FILE__

cdef extern from "Python.h":
    int __LINE__




"""
#In Cython we would do something equivalent to CoolProp.pxd 
# This class inheritate from State_Flooded which calculates the mixture properties 
# and returns the values in a similar way of State class in CoolProp.pxd
"""

#TODO:  check consistency with StateFlooded.pyx   get_h_m or get_h. Maybe we should differentiate from standard State class and use _m to indicate mixture properties

cdef class State_Flooded:
    cdef AbstractState pAS
    cdef readonly bytes Fluid, phase
    cdef int iFluid,iParam1,iParam2,iOutput
    cdef int iLiquid    # Maybe ??
    
    cdef double T_, rho_, p_, xL_,
    
    cpdef set_Fluid(self, string Fluid, string backend)
    cpdef set_Liquid(self, string Fluid)  #No need for backend
    
    #Now the properties are all defined above and we can be consistend with the definition
    #I used in containers_20150419  @properties
    
    cpdef speed_test(self, int N)
    cpdef update(self, dict params)
    cpdef update_ph(self, double p, double h)
    cpdef update_Trho(self, double T, double rho)
    cpdef State copy(self)
    cpdef double Props(self, constants_header.parameters iOutput) except *
    cpdef long Phase(self) except *
    cpdef double get_Q(self) except *
    cpdef double get_T(self) except *
    cpdef double get_p(self) except *
    cpdef double get_h_m(self) except *
    cpdef double get_rho_m(self) except *
    cpdef double get_s_m(self) except *
    cpdef double get_u_m(self) except *
    cpdef double get_visc_m(self) except *
    cpdef double get_cond_m(self) except *
    cpdef double get_cp_m(self) except *
    cpdef double get_cp0(self) except *   # TO CHECK
    cpdef double get_cv_m(self) except *
    cpdef double get_dpdT(self) except *
