#cython: cdivision=True
    
from libc.math cimport cos, sin, sqrt, M_PI as pi

from PDSim.misc._listmath import listm
from PDSim.misc._listmath cimport listm

from CoolProp.State import State
from CoolProp.State cimport State

cdef class _Recip(object):
    
    cpdef dict __cdict__(self):
        return dict(crank_length = self.crank_length, 
                    connecting_rod_length = self.connecting_rod_length, 
                    A_piston = self.A_piston, 
                    omega = self.omega, 
                    piston_diameter = self.piston_diameter,
                    V_dead = self.V_dead)
        
    cpdef tuple V_dV(self, double theta):
        cdef double x_calc, x_2, V, dV
        x_calc = self.crank_length*cos(theta) + sqrt(self.connecting_rod_length**2 - self.crank_length**2*sin(theta)**2)
        x_2 = (self.connecting_rod_length + self.crank_length)-x_calc
        V=x_2*self.A_piston+self.V_dead
        dV=-(-self.crank_length*sin(theta) - (self.crank_length**2*sin(2*theta))/
             (2*sqrt(self.connecting_rod_length**2 - self.crank_length**2*sin(theta)**2)))*self.A_piston
        return V, dV
    
    cpdef listm heat_transfer_callback(self, double theta):
        cdef double T_w, V,dV,D_h,A_ht,Pr,rho,k,mu,T,u,Re,h_c,Q,cp
        cdef State S
        
        T_w  = self.Tlumps[0] #[K]
        V,dV = self.V_dV(theta) #[m3,m3/radian]
        D_h = self.piston_diameter #[m]
        A_ht = pi*self.piston_diameter*(V/self.A_piston) #[m2]
        
        S = self.CVs['A'].State
        cp = S.get_cp() #[kJ/kg/K]
        rho = S.get_rho() #[kg/m3]
        k = S.get_cond() #[kW/m-K]
        mu = S.get_visc() #[Pa-s]
        T = S.get_T() #[K]
        Pr = cp*mu/k
        u = abs(0.5*dV*self.omega/self.A_piston) #[m/s]
        Re = (rho*u*D_h)/mu #[kg/m3*m/s*m/Pa/s]=[-]
        
        h_c = 0.053*(k/D_h)*Pr**(0.6)*Re**(0.8) #[kW/m2/K]
        Q = h_c*A_ht*(T_w-T)   #Net heat into control volume [kW]
        if self.CVs.N > 1:
            return listm([Q,0])
        else:
            return listm([Q])