""""
Davide Ziviani 2/16/2016

RollingPiston
"""

##########################################
###          Part 1. Imports           ###
##########################################

#Here we import some python packages
from __future__ import division
from math import pi, cos, sin, asin,tan,sqrt
from scipy.constants import g
from scipy.integrate import quad,simps
from time import clock
import numpy as np
import pylab
import os, sys

# If the following line is uncommented, python will try to use a local version
# of PDSim.  This is handy for debugging purposes.  Generally you want this line 
# commented out
# PDSim should also be built using a command like python build_ext --inplace to keep all the extension modules next to the .pyx files
#sys.path.insert(0, os.path.abspath('..'))

#Here we import the things from PDSim we need
from PDSim.flow.flow import FlowPath
from PDSim.flow import flow_models
from PDSim.misc.datatypes import arraym,listm
from PDSim.core.containers import ControlVolume, Tube
from PDSim.core.core import PDSimCore
from PDSim.plot.plots import debug_plots

#Imports from CoolProp (fluid property database)
from CoolProp import State,AbstractState
from CoolProp import CoolProp as CP
    
############################################################
###   Part 2. Declaration of Rolling Piston Compressor   ###
###########################################################

class RollingPistonCompressor(PDSimCore):
    
    def __init__(self):
        #Initialize the base class that PURecip is derived from
        PDSimCore.__init__(self)
        

    def V_dV_ParkYC(self,theta):
        
        """
        Park Y.C. 'Transient analysis of a variable speed rotary compressor',Energy Conversion and Management 2010 
        
        Compression Chamber
        
        """
        a = self.Rr/self.Rc
        V = self.Hc/2*self.Rc**2*((1-a**2)*theta - (1-a)**2*(sin(2*theta))/2 - a**2*asin( (1/a-1)*sin(theta)) -a*(1-a)*sin(theta)*sqrt(1-(a-1)**2*(sin(theta))**2)) - self.b/2*self.Hc*self.Rc*( 1-(1-a)*cos(theta)- sqrt((1-a)**2*(cos(theta))**2+2*a-1) )
        
        dVdtheta = (self.Hc*self.Rc**2/2)* (-((1/a-1)*a**2*cos(theta))/(sqrt(1-(1/a-1)**2*(sin(theta))**2)) - a**2 + (1-a)**2*(-cos(2*theta))-a*(1-a)*cos(theta)*sqrt(1-(1/a-1)**2*(sin(theta))**2)+ ((1/a-1)**2*a*(1-a)*(sin(theta))**2*cos(theta)  )/(sqrt(1-(1/a-1)**2*(sin(theta))**2)) + 1   ) - (self.b*self.Hc*self.Rc/2)*( (1-a)*sin(theta)+ ((1-a)**2*sin(theta)*cos(theta))/(sqrt((1-a)**2*(cos(theta))**2+2*a-1))  )
        
        return abs(V),dVdtheta


    def Vc_dVc_ParkYC(self,theta):
    
        V,dVdtheta = self.V_dV_ParkYC(2*pi - theta)
        
        if V < 0.0:
            raise Exception('V is negative!')

        return V,dVdtheta
    
    def x_theta(self,x):
        """
        Vane displacement with the crankshaft angle
        """
        a = self.Rr/self.Rc
        return self.Rc*(1-(1-a)*cos(x) - sqrt((1-a)**2*(cos(x))**2+2*a-1 ) )
    
    def V_disp(self,theta_s,theta_d):
        """
        Zheng et al. 'Experimental verification of a rolling-piston expander that
        applied for low-temperature ORC', Applied Energy 112(2013)1265-1274.
        
        ParkY.C. 'Transient analysis of a variable speed rotary compressor' 
        Energy Conversion and Management 51(2010) 277-287 
        
        """
        V_theta_s = self.V_dV_ParkYC(theta_s)[0]
        V_theta_d = self.Vc_dVc_ParkYC(2*pi - theta_d)[0]
        Int_x_dtheta, err = quad(self.x_theta,0,2*pi)
        Vdisp = pi*self.Hc*(self.Rc**2 - self.Rr**2) - self.b*self.Hc*Int_x_dtheta - (V_theta_s + V_theta_d)
        
        return Vdisp
    
    
    def A_HT(self,theta):
        """
        Calculation of chamber total heat transfer area
        """
        a = self.Rr/self.Rc
        
        hv_theta = self.Rc*(1-(1-a)*cos(theta)-sqrt((1-a)**2*(cos(theta))**2) +2*a-1)
        f_theta = (1-a**2)*theta-0.5*(1-a)**2*sin(2*theta) - a**2*asin((1/a-1)*sin(theta)) - a*(1-a)*sin(theta)*sqrt(1-(1/a-1)**2*(sin(theta))**2)
        
        phi = pi/2 + theta - asin((self.Rc - hv_theta + self.e*cos(pi - theta))/self.Rr )
        
        #delta = self.Rc- (2*self.e*cos(theta-phi)+sqrt(4*e**2*(cos(theta-phi))**2-4*(e**2 - self.Rr**2)))/2
        
        A = theta*(self.Rc + self.Rr)*self.Hc+ 2*0.5*self.Rc**2*f_theta
        return A
    
    def A_HT_c(self,theta):
    
        A = self.A_HT(2*pi - theta)
        return A  
    
    def r_chamber_mean(self,theta):
        """
        Calculates the mean chamber radius to be used for heattransfer coefficient calculation
        """
        a = self.Rr/self.Rc
        N = 100
        phi_array = np.linspace(0,theta,N)    
        R_array = np.zeros(len(phi_array))
        ravg_array = np.zeros(len(phi_array))
        
        for i in range(len(phi_array)):
        
            R_array[i] =  self.Rc*((1-a)*cos(theta-phi_array[i])+sqrt((1-a)**2*(cos(theta- phi_array[i]))**2+(2*a-1)))
            ravg_array[i] = (self.Rc + R_array[i])/2
        
        r = np.sum(ravg_array)/N
        
        return r
        
        
    def Suction(self,FlowPath):
        
        if FlowPath.key_up=='A':
            # pressure in compressor higher than the inlet line
            # valve is closed - no flow

            return 0.0
        else:
            try:
                FlowPath.A=self.A_suction
                mdot=flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
                return mdot
            except ZeroDivisionError:
                return 0.0

        
    def Discharge(self,FlowPath):
            if FlowPath.key_down=='B':
                ## pressure in compressor lower than the discharge line
                ## valve is closed - no flow
                return 0.0
            else:
                try:
                    FlowPath.A=self.A_discharge
                    mdot=flow_models.IsentropicNozzle(FlowPath.A,
                                                    FlowPath.State_up,
                                                    FlowPath.State_down)
                    return mdot
                except ZeroDivisionError:
                    return 0.0
            
    def TubeCode(self, Tube):
        """ 
        A thin wrapper of the isothermal wall tube from flow_models.py 
        """
        Tube.Q = flow_models.IsothermalWallTube(Tube.mdot, 
                                                Tube.State1, 
                                                Tube.State2,
                                                Tube.fixed, 
                                                Tube.L, 
                                                Tube.ID,
                                                T_wall = self.Tlumps[0])
    
    def heat_transfer_callback(self, theta):
        """
        A callback used by PDSimCore.derivs to calculate the heat transfer
        to the gas in the working chamber.
        We return an arraym instance the same length as the number of CV in existence

        """
        if self.HT_on == True:
            Qarray = []

            for key in self.CVs.exists_keys:
                if key == 'A':
                    #Suction chamber
                    V,dVdtheta = self.V_dV_ParkYC(theta)
                    A_ht = self.A_HT(theta) #[m2]
                    
                elif key == 'B':
                    #Compression chamber
                    V,dVdtheta = self.Vc_dVc_ParkYC(theta)
                    A_ht = self.A_HT_c(theta) #[m2]
                else:
                    raise Exception('Something wrong with CVs names')
                
                #Avoid division by zero when Aht is zero
                try:    
                    T_w  = self.Tlumps[0] #[K]
                    r_avg = self.r_chamber_mean(theta)
                    
                    D_h = 4*V/A_ht #[m]
                    
                    Pr = self.CVs[key].State.Prandtl #[-]
                    rho = self.CVs[key].State.rho #[kg/m3]
                    k = self.CVs[key].State.k #[kW/m-K]
                    mu = self.CVs[key].State.visc #[Pa-s]
                    T = self.CVs[key].State.T #[K]
                    u = abs(0.5*self.omega*self.Rr) #[m/s]
                    Re = (rho*u*D_h)/mu #[kg/m3*m/s*m/Pa/s]=[-]
                    h_c = 0.023*(k/D_h)*Pr**(0.4)*Re**(0.8)*(1.0+1.77*(D_h/r_avg)) #[kW/m2/K]
                    Q = h_c*A_ht*(T_w-T)   #Net heat into control volume [kW]
                except:
                    Q = 0.0
                    
                Qarray.append(Q)

            return arraym(Qarray)
        
        elif self.HT_on == False:
            return arraym([0.0]*self.CVs.N)
            
        else:
            raise Exception('HT_on is either True or False')

        
    def mechanical_losses(self):
        """
        To be added
        """
        
        return self.Wdot_parasitic


    def friction_losses(self):
        #TODO: to be complete
        """
        Friction losses are caused by rubbing parts
        
        Ooi K.T. 'Design optimization of a rolling piston compressor 
        for refrigeration', Applied Therm. Eng. 25(2005) 813-829
    
        """
        DELTA_omega = self.omega - omega_r
        #Frictional loss between eccentric and the inner surface of the roller Ler
        
        Ler = DELTA_omega*(2*pi*eta*DELTA_omega)*self.Hc*R_e**3/delta2
        
        #Frictional loss between the roller face and the cylinder head faces Lrc
        Lrc = self.omega*(2*pi*omega_r*eta*(Rr**4-Re**4))/delta_1r
        
        #Frictional loss between eccentric face and the cylinder head faces Lec
        Lec = self.omega*(pi*self.omega*eta*(2*Re**4-Rs**4))/delta3
        
        #Frictional loss between vane tip and roller Lv
        
        Lv = Fvt*Rr*omega_r + (e*self.omega*cos(theta))/cos(alpha)
        
        #Frictional loss between vane sides and vane slot Ls
        Ls = (Ft1 + Ft2)*xdot
        
        #Frictional loss between shaft and bearing S, Lbs
        Lbs = self.omega*(2*pi*eta*self.omega*lbs*Rs**3)/delta4
        
        #Frictional loss between shaft and beating L, Lbl
        Lbl = self.omega*(2*pi*eta*self.omega*lbl*Rs**3)/delta4
        
        #self.Wdot_friction = Ler + Lrc + Lec
        self.Wdot_friction = 0.0
        
        return self.Wdot_friction
    
    def ambient_heat_transfer(self):
        """
        The ambient heat transfer for the compressor in kW
        
        Returns a positive value if heat is added to the compressor from the 
        ambient
        
        The heat transfer is obtained for natural convection from a vertical cylinder
        """
        ref = 'Air'
        pamb = 101.325 #[kPa]
        Tfilm = (self.Tamb + self.Tlumps[0])/2 #[C]
        HEOS = AbstractState("HEOS",ref)
        HEOS.update(CP.PT_INPUTS,pamb*1000,Tfilm +273.15)
        rho_film = HEOS.rhomass() #[kg/m3]
        mu_film = HEOS.viscosity() #[Pa-s]
        nu_film = mu_film/rho_film #[m2/s]
        k_film = HEOS.conductivity()/1000 #[kW/m2-K]
        Pr_film = HEOS.Prandtl() #[-]
        cp_film = HEOS.cpmass()/1000 #[kJ/kg-K]
        beta_film = 1 #Thermal expansion 
        Gr_L = g*beta_film*(Tfilm - self.Tamb)*self.Hshell**3/nu_film**2
        Ra_L = Gr_L*Pr_film
        
        Nu_0 = 0.68 #for vertical plate
        h_shell = (k_film/self.Hshell)*(Nu_0**(1/2)+0.387*Ra_L**(1/6)*(1+(0.492/Pr_film)**(9/16))**(-8/27)   )**2 #kW/m2/K
        return h_shell*self.A_shell*(self.Tamb-self.Tlumps[0]) #[kW]
        
    def lump_energy_balance_callback(self):
        """
        A callback used in PDSimCore.solve to do the energy balance on the lump
        
        Note: we neglect heat transfer to the gas in the working chamber
        """
        #Mechanical losses are added to the lump
        self.Wdot_mechanical = self.mechanical_losses() #[kW]
        #Friction losses
        self.Wdot_friction = friction_losses(self) #[kW]
        #Heat transfer between the shell and the ambient
        self.Qamb = self.ambient_heat_transfer() #[kW]
        return self.Wdot_mechanical + self.Qamb

    def step_callback(self,theta,h,Itheta):
        self.theta = theta
        return False, h
 
######################################################
###    Part 3. Execution of RollingPistonCompressor###
######################################################
        
def GeometryPlot():
    rolling=RollingPistonCompressor()                               #Instantiate the model
    rolling.e = 5.0e-3                                                #Eccentricity [m]
    rolling.Rc = 27.0e-3                                              #Radius of cylinder [m3]
    rolling.Rr = 22.0e-3                                              #Radius roller [m]
    rolling.rv = 2.0e-3                                               #Radius of vane tip [m]
    rolling.Hc =  24.0e-3                                             #Compressor cylinder height [m]
    rolling.b =  5.0e-3                                               #Vane thickness [m]
    rolling.hv = 24.0e-3                                               #Vane height [m]

    
    N = 100
    theta = np.linspace(0,2*pi,N,endpoint = True)
    V = np.zeros(len(theta))
    V2 = np.zeros(len(theta))
    dV = np.zeros(len(theta))
    dV2 = np.zeros(len(theta))
    V_Cai = np.zeros(len(theta))
    V_ParkYC = np.zeros(len(theta))
    dV_ParkYC = np.zeros(len(theta))
    Vs_ParkYC = np.zeros(len(theta))
    dVs_ParkYC = np.zeros(len(theta))
    A_ht = np.zeros(len(theta))
    A_ht_c = np.zeros(len(theta))
    
    for i in range(len(theta)):
        V_ParkYC[i],dV_ParkYC[i]=rolling.Vc_dVc_ParkYC(theta[i])
        Vs_ParkYC[i],dVs_ParkYC[i]=rolling.V_dV_ParkYC(theta[i])


    f, (ax1) = pylab.subplots(1, 2,figsize=(10,6))
    f.subplots_adjust(bottom=0.15,left= 0.1)  
    ax1.plot(theta,V_ParkYC*1e6,'r--',lw = 2,label = 'Vc Park Y.C.')
    ax1.plot(theta,dV_ParkYC*1e6,'r-.',lw = 2,label = 'dVc Park Y.C.')
    ax1.plot(theta,Vs_ParkYC*1e6,'g--',lw = 2,label = 'Vs Park Y.C.')
    ax1.plot(theta,dVs_ParkYC*1e6,'g-.',lw = 2,label = 'dVs Park Y.C.')
    ax1.legend(loc = 'best', fontsize = 10)
    ax1.set_xlabel(r'$\theta [rad]$',fontsize = 15)
    ax1.set_ylabel(r'$V \, [cm^3]$',fontsize = 15)
    ax1.set_title(r'V($\theta$)', fontsize = 15)
    ax1.set_ylim(0,20)

    pylab.show()
    

def RollingCompressor():

    rolling=RollingPistonCompressor()                               #Instantiate the model
    rolling.e = 5.0e-3                                                #Eccentricity [m]
    rolling.Rc = 27.0e-3                                              #Radius of cylinder [m3]
    rolling.Rr = 17.0e-3                                              #Radius roller [m]
    rolling.rv = 2.0e-3                                               #Radius of vane tip [m]
    rolling.Hc =  24.0e-3                                             #Compressor cylinder height [m]
    rolling.b =  5.0e-3                                               #Vane thickness [m]
    rolling.h = 24.0e-3                                               #Vane height [m]
    rolling.Hshell = 30e-3                                          #Height of compressor shell
    rolling.Rshell = 35e-3                                          #compressor shell outer diameter
    rolling.Nmot = 3000                                                #Electric motor rotational speed
    rolling.omega = 2*pi*rolling.Nmot/60                                          #Frequency, rad/sec (60Hz)
    rolling.d_discharge=0.005;                                     #discharge port diameter [m]
    rolling.d_suction=0.01;                          #suction diameter [m]
    rolling.A_discharge=pi*rolling.d_discharge**2/4
    rolling.A_suction=pi*rolling.d_suction**2/4
    #Suction angles
    rolling.thetad = 10*pi/180
    rolling.beta = 50*pi/180
    #Discharge angles
    rolling.thetas = 28*pi/180
    rolling.psi = 20*pi/180
    
    #Selection of HT model
    rolling.HT_on = True    
    #These are parameters needed for the ambient heat transfer model
    rolling.A_shell = pi*rolling.Hshell*2*rolling.Rshell #[m2]
    rolling.Tamb = 25 + 273.15                    #[K] 
    
    #Parameters for the mechanical losses model (simplified)
    rolling.Wdot_parasitic = 0.01         #Parasitic losses [kW]
    
    Ref='R410A'
    pin = 1087
    Tin = 33.2+273.15
    pc = 2500
    
    inletState=State.State(Ref,dict(T=Tin, P = pin))
    T2s = rolling.guess_outlet_temp(inletState,pc)

    outletState=State.State(Ref,{'T':T2s,'P':pc})
    
    #Guess mass flow rate
    rolling.Vdisp = rolling.V_disp(rolling.thetas,rolling.thetad) #Displacement volume
    mdot_guess = inletState.rho*rolling.Vdisp*rolling.Nmot/60
    
    #print 'mdot_guess [kg/s]:',mdot_guess
    #print 'Vdisp [cm3]:',rolling.Vdisp*1e6
    
    
    
    #Add the inlet tube
    rolling.add_tube( Tube(key1='inlet.1',
                        key2='inlet.2',
                        L=0.03,ID=0.02,
                        mdot=mdot_guess, 
                        State1=inletState.copy(),
                        fixed=1,
                        TubeFcn=rolling.TubeCode) 
                        )
    
    #Add the outlet tube
    rolling.add_tube( Tube(key1='outlet.1',
                            key2='outlet.2',
                            L=0.03,ID=0.02,
                            mdot=mdot_guess, 
                            State2=outletState.copy(),
                            fixed=2,
                            TubeFcn=rolling.TubeCode) 
                            )
    
    #Add the control volumes.
    """
    'A' = Suction CV
    'B' = Compression/discharge CV
    """
    rolling.add_CV( ControlVolume(key='A',
                               initialState=inletState.copy(),
                               VdVFcn=rolling.V_dV_ParkYC,
                               becomes='A') )
    rolling.add_CV( ControlVolume(key='B',
                               initialState=outletState.copy(),
                               VdVFcn=rolling.Vc_dVc_ParkYC,
                               becomes='B') )

    
    #Add the flow paths that link flow nodes together
    rolling.add_flow(FlowPath(key1='inlet.2',key2='A',MdotFcn=rolling.Suction))
    rolling.add_flow(FlowPath(key1='outlet.1',key2='B',MdotFcn=rolling.Discharge))
    
    rolling.connect_callbacks(step_callback = rolling.step_callback,
                            endcycle_callback=rolling.endcycle_callback, # Provided by PDSimCore
                            heat_transfer_callback=rolling.heat_transfer_callback,
                            lumps_energy_balance_callback = rolling.lump_energy_balance_callback
                            )
    
    t1=clock()
    rolling.solve(key_inlet='inlet.1',
                key_outlet='outlet.2',
                solver_method = 'Euler',
                OneCycle = False,
                UseNR = False,
                #cycle_integrator_options = dict(tmax = 2*pi)
                )
    print 'time taken',clock()-t1,'s'
    
    debug_plots(rolling)
    
if __name__=='__main__':    
    
    GeometryPlot()
    RollingCompressor()
    