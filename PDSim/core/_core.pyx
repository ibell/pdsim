from __future__ import division

import cython
cimport cython

#Import the listm type used in PDSim
from PDSim.misc._listmath import listm
from PDSim.misc._listmath cimport listm

cdef class _PDSimCore:
        
    cpdef _derivs(self,double theta, listm x, heat_transfer_callback=None):
        cdef listm xL,V,dV,rho,dpdT,h,cv,p,m,summerdxL,summerdT,dmdtheta
        #1. Calculate the volume and derivative of volume of each control volumes
        #    Return two lists, one for the volumes, second for the derivatives of volumes 
        V,dV=self.CVs.volumes(theta)
    
        if self.__hasLiquid__==True:
            raise NotImplementedError('Not coded')
        else:
            self.CVs.updateStates('T',x[0:self.CVs.Nexist],'D',x[self.CVs.Nexist:2*self.CVs.Nexist])
            
        #2. Calculate the mass flow terms between the control volumes
        self.Flows.calculate(self)
        summerdT,summerdm=self.Flows.sumterms(self)
        
        #3. Calculate the heat transfer terms
        if heat_transfer_callback is not None:
            Q=heat_transfer_callback(theta)
            if not len(Q) == self.CVs.Nexist:
                raise ValueError
        else:
            Q=0.0
        
        ## Calculate properties and property derivatives
        ## needed for differential equations
        rho = listm(self.CVs.rho)
        v = 1.0/rho
        m = V*rho
        T = listm(self.CVs.T)
        h = listm(self.CVs.h)
        cv = listm(self.CVs.cv)
        dpdT = listm(self.CVs.dpdT)
        p = listm(self.CVs.p)
        
        self.V_=V
        self.dV_=dV
        self.rho_=rho
        self.m_=m
        self.p_=p
        self.T_=T
                
        dudxL=0.0
        summerdxL=0.0*T
        xL=0.0*T
        
        dmdtheta=summerdm
        dxLdtheta=1.0/m*(summerdxL-xL*dmdtheta)
        dTdtheta=1/(m*cv)*(-1.0*T*dpdT*(dV-v*dmdtheta)-m*dudxL*dxLdtheta-h*dmdtheta+Q/self.omega+summerdT)
        drhodtheta = 1.0/V*(dmdtheta-rho*dV)
        
        if self.__hasLiquid__==True:
            raise NotImplementedError('Not Coded')
        else:
            return listm(list(dTdtheta)+list(drhodtheta))