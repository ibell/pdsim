import cython
cimport cython

#Import the listm type used in PDSim
from PDSim.misc._listmath cimport listm
from PDSim.flow._Flow cimport _FlowPath

cdef class _PDSimCore:
    cdef bint __hasLiquid__
    cdef readonly int NCV
    cdef public listm V_,dV_,T_,p_,rho_,summerdm_,summerdT_,summerdxL_
    cdef list exists
    
    cpdef _derivs(self,double theta,listm x,heat_transfer_callback=*)
    #@cython.locals(V=np.ndarray,dV=np.ndarray,I=cython.int,icv=cython.int,T=np.ndarray,rho=np.ndarray,xL=np.ndarray)
    #cpdef derivs(self, double theta, x, heat_transfer_callback=*)
##     
##     #@cython.locals(disableAdaptive=cython.bint,t0=cython.double,tmax=cython.double)
    #cpdef cycle(self, double hmin=?,double tmin=?,double tmax=?,double eps_allowed=?,double step_relax=?,step_callback=?,heat_transfer_callback=?)