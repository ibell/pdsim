"""
From Flores 2006, JOURNAL BEARINGS SUBJECTED TO DYNAMIC LOADS: THE ANALYTICAL MOBILITY METHOD
"""

from __future__ import division
import numpy as np
from math import pi
import matplotlib.pyplot as plt
# If scipy is available, use its spline interpolation function, otherwise, 
# use our implementation (for packaging purposes)
try:
    import scipy.interpolate as interp
except ImportError:
    import PDSim.misc.scipylike as interp

if __name__=='__main__':
    omega_bar = 2000/60.0*2*pi
    theta = np.array([-0.570621, 15.768046, 34.544827, 51.431137, 68.188551, 87.101887, 104.527875, 120.872984, 139.924246, 156.677399, 173.557215, 192.464921, 209.346919, 226.090992, 244.595321, 261.887695, 278.912488, 297.704183, 314.605203, 331.368349, 350.292845, 367.055839, 386.113492, 403.005533, 419.900923, 438.686784, 455.442879, 472.470614, 491.525375, 508.148313, 525.18665, 544.253434, 560.621675, 577.93155, 596.452872, 613.759501, 630.659151, 649.570865, 666.186194, 683.213269, 701.726729, 718.618619])/180.0*pi
    t = theta/omega_bar
    F = np.array([-2.012326, -1.595861, -1.365927, -1.182637, -1.204212, -1.020913, -0.809306, -0.604367, -0.512673, -0.394341, 0.00214, 0.370316, 0.695178, 1.111646, 1.413198, 1.574837, 1.6482, 1.38846, 1.088737, 0.878955, 0.695831, 0.491044, 0.372878, 0.367959, 0.253114, 0.184913, 0.206643, 0.183404, 0.160174, 0.116946, -0.254394, -0.672362, -1.226917, -1.639896, -1.896306, -2.202689, -2.457441, -2.220844, -2.014238, -2.015825, -2.014074, -2.013995])*1000
    
    F_interp = interp.splrep(t, F, k = 2, s = 0)
    
    mu = 0.00416 #[Pa-s]
    L = 0.0254 #[m]
    D = 0.0635 #[m]
    c = 35.56e-6 #[m]
    L_over_D = L/D #[-]
    c_R = c/(D/2.0) #[-]
    
    from PDSim.core.bearings import journal_bearing
    
    def f(x,F_t,t):
        
        epsilon = x[0]
        phi = x[1]
            
        zeta = epsilon*np.cos(phi)
        kappa = epsilon*np.sin(phi)
        
        if F_t > 0:
            M_zeta = (1-zeta)**2.5/(pi*L_over_D**2)
            M_kappa = -4*kappa*(1-zeta)**1.5/(pi**2*L_over_D**2)
        else:
            M_zeta = (1+zeta)**2.5/(pi*L_over_D**2)
            M_kappa = 4*kappa*(1+zeta)**1.5/(pi**2*L_over_D**2)
        
        M_epsilon = M_zeta*np.cos(phi)+M_kappa*np.sin(phi)
        M_phi = -M_zeta*np.sin(phi)+M_kappa*np.cos(phi)
        
        epsilon_dot = F_t*(c_R)**2/(mu*L*D)*M_epsilon
        phi_dot = F_t*(c_R)**2/(mu*L*D*epsilon)*M_phi+omega_bar
        
        Fb = F_t*c*epsilon/D*np.sin(phi)+2*pi*mu*omega_bar*L/(c*np.sqrt(1-epsilon**2))*(D/2)**2
        f = Fb/F_t
        
        print(f)
        
        #plt.plot(t,c*(1-epsilon)*1e6,'o')
    #    plt.plot(kappa,-zeta,'o')
    #    plt.plot(t,F_t,'o')
        plt.plot(t,f,'o')
        #plt.plot(t,phi,'o')
        
    #    plt.plot(t,c*(1-epsilon)*1e6,'o',t,journal_bearing(r_b = D/2,L  = L, eta_0 = mu, omega = omega_bar, c = c, W = abs(F_t))['h_min']*1e6,'.')
        
        
        return np.array([epsilon_dot, phi_dot])
    
    x = np.array([0.1, 1.5])
    t = 0
    dt = 0.00005
    
    while t < 0.12:
        
        # Interpolate to find F
        F_t = interp.splev(t%(4*pi/omega_bar), F_interp)
        
        # Get the derivatives
        xdot = f(x, F_t, t)
        
        # Euler step
        x += dt*xdot
        
        # Update t
        t += dt
        
    #plt.xlim(-1,1)
    #plt.ylim(-1,1)
    #plt.gca().set_aspect(1)
    plt.show()