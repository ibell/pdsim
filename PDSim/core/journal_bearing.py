from __future__ import division,print_function

import numpy as np
from math import pi, atan
# If scipy is available, use its optimization function, otherwise, 
# use our implementation (for packaging purposes)
try:
    import scipy.optimize as optimize
except ImportError:
    import PDSim.misc.scipylike as optimize
    
import matplotlib.pyplot as plt

N = 61
#e_mat=[0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9];

phi_star = pi

def TwoDGriddedIntegrate(I,N):
    # Average the center of each cell based on its neighboring nodes
    return np.sum(np.sum((I[0:N-1,0:N-1]+I[1:N,0:N-1]+I[0:N-1,1:N]+I[1:N,1:N])))/4
    
def TwoDGriddedIntegrate2(PHI,Y,I):
    
    #Integrate along phi direction for each y, then do a trapezoidal integration of each of the y
    plt.plot(Y[1,:],np.trapz(I,PHI,axis = 0))
    plt.show()
    return np.trapz(np.trapz(I,PHI,axis = 0),Y[1,:])

def OBJECTIVE(phi_star, epsilon, plot = False, output = False):

    PHI = np.tile(np.linspace(0,phi_star,N).T,(N,1)).T
    Y = np.tile(np.linspace(0,1,N),(N,1))
    
    dPHI = phi_star/(N-1)
    dY = 1/(N-1)
    sinPHI=np.sin(PHI)
    P = 0*PHI
    Pnew = 0*PHI
    f = 0*PHI
    df = 0*PHI
    _lambda = 1
    change = 999    
    
    eps=1e-6;
    count=0;
    while (change>eps):
        
        #Calculate geometric parameters
        H=1+epsilon*np.cos(PHI);
        H3=H**3;
        
        #Coefficients
        A=H3[2:N,1:N-1]
        B=H3[0:N-2,1:N-1]
        C=H3[1:N-1,1:N-1]

        #Calculate residuals
        f[1:N-1,1:N-1] = -(4*A+4*B+2*_lambda*dPHI**2/dY**2*C)*P[1:N-1,1:N-1]+(3*A+B)*P[2:N,1:N-1]+(A+3*B)*P[0:N-2,1:N-1]+(_lambda**2*dPHI**2/dY**2*C)*(P[1:N-1,2:N]+P[1:N-1,0:N-2])+24*dPHI**2*epsilon*sinPHI[1:N-1,1:N-1]

        #Calculate derivative
        df[1:N-1,1:N-1]=-(4*A+4*B+2*_lambda*dPHI**2/dY**2*C);

        #Evaluate P_new=P_old-f/dfdP
        P[1:N-1,1:N-1]=P[1:N-1,1:N-1]-f[1:N-1,1:N-1]/df[1:N-1,1:N-1];

        #Evaluate change
        change=np.max(np.max(np.abs(f[1:N-1,1:N-1]/df[1:N-1,1:N-1])));

        if count % 1000 == 0:
            print(change)
            
        count += 1
        
    if output:
        Wx=dY*dPHI*np.sum(np.sum(np.sin(PHI)*P))
        Wz=-dY*dPHI*np.sum(np.sum(np.cos(PHI)*P))
        Wr = np.sqrt(Wx**2+Wz**2)
        PHI_angle = atan(Wx/Wz)
        B_j = 1/(pi*Wr)
        DPDPHI = 0*Y
        DPDPHI[0:N-2,0:N] = (P[1:N-1,0:N]-P[0:N-2,0:N])/(dPHI)
        DPDPHI[N-1:N-1,0:N] = (P[N-1:N,0:N]-P[N-2:N-2,0:N])/(dPHI)
        
        integrand = 1/H
        #integrand = H/2*DPDPHI+1/H
        
        Fb1 = dPHI*dY*np.sum(np.sum(integrand))
        Fb2 = dPHI*dY*TwoDGriddedIntegrate(integrand,N)
        Fb3 = TwoDGriddedIntegrate2(PHI,Y,integrand)
        mu_rb_c = Fb3/Wr # mu*r_b/c
        print('Fb1,Fb2,Fb3',Fb1,Fb2,Fb3)
        print('B_j', B_j)
        print('mu*rb/c', mu_rb_c)
        #print 'mu*rb/c', mu_rb_c/12.8
        print('PHI_angle', PHI_angle/pi*180)
        plt.contour(PHI,Y,H/2*DPDPHI+1/H)
        plt.show()
    
    if plot:
        plt.contour(PHI,Y,P,30)
        plt.show()
    
    return np.sum(3*P[N-1,N//2+1]-4*P[N-2,N//2+1]+P[N-3,N//2+1])/(2*dPHI)
        
if __name__=='__main__':
    print(optimize.newton.__doc__); quit()
    phi_star = optimize.newton(OBJECTIVE, pi, args = (0.6,), tol = 0.004)
    
    OBJECTIVE(phi_star,0.6,plot = True, output = True)
