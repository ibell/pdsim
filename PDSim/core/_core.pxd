import cython
cimport cython

#Import the listm type used in PDSim
from PDSim.misc._listmath cimport listm

cdef class _PDSimCore:
    pass
    #cdef bint __hasLiquid__
    #cdef readonly int NCV
    #cdef public listm V_,dV_,T_,p_,rho_,summerdm_,summerdT_,summerdxL_
    #cdef list exists