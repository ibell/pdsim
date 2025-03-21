from __future__ import division, print_function, absolute_import
from math import pi, cos, sin, asin,tan,sqrt,log,log10
from scipy.constants import g
from scipy.integrate import quad,simps
import numpy as np
import os, sys

from PDSim.flow.flow import FlowPath
from PDSim.flow import flow_models
from PDSim.misc.datatypes import arraym,listm
from PDSim.core.containers import ControlVolume, Tube
from PDSim.core.core import PDSimCore
from PDSim.heattransfer.HTC import HTC
from PDSim.plot.plots import debug_plots
from PDSim.flow.flow_models import ValveModel
from PDSim.rolling import rolling_geo

#Imports from CoolProp (fluid property database)
from CoolProp import State,AbstractState
from CoolProp import CoolProp as CP
import CoolProp
  
class struct(object):
    pass

class RollingPiston(PDSimCore):
    
    def __init__(self):
        #Initialize the base class that RollingPiston is derived from
        PDSimCore.__init__(self)
        
        ## Define the geometry structure
        self.geo=rolling_geo.geoVals()

    def set_rolling_geo(self,e,Rc,Rc_outer,Rr,Rs,rv,Hc,b,hv,offsetRoller,D_sp,D_dpc,theta_sm,theta_dm,delta_theta_s,delta_theta_d):

        self.geo.e = e     
        self.geo.Rc = Rc  
        self.geo.Rc_outer = Rc_outer                   
        self.geo.Rr = Rr
        self.geo.Rs = Rs    
        self.geo.rv = rv    
        self.geo.Hc = Hc   
        self.geo.b =  b    
        self.geo.hv = hv 
        self.geo.offsetRoller = offsetRoller  
        self.geo.theta_sm = theta_sm
        self.geo.theta_dm = theta_dm
        self.geo.D_sp  = D_sp
        self.geo.D_dpc = D_dpc
        self.geo.theta_sm = theta_sm
        self.geo.theta_dm = theta_dm
        self.geo.delta_theta_s = delta_theta_s
        self.geo.delta_theta_d = delta_theta_d
        
        #Set the flags to ensure all parameters are fresh
        self.__Setrolling_geo__=True

    def V_dV_ParkYC(self,theta):
        
        """
        Park Y.C. 'Transient analysis of a variable speed rotary compressor',Energy Conversion and Management 2010 
        Mathematical description of compression chamber
        """
        a = self.Rr/self.Rc
        Vclearance = 1e-7 #m^3
        V = Vclearance +self.Hc/2*self.Rc**2*((1-a**2)*theta - (1-a)**2*(sin(2*theta))/2 - a**2*asin( (1/a-1)*sin(theta)) -a*(1-a)*sin(theta)*sqrt(1-(a-1)**2*(sin(theta))**2)) - self.b/2*self.Hc*self.Rc*( 1-(1-a)*cos(theta)- sqrt((1-a)**2*(cos(theta))**2+2*a-1) )
        
        dVdtheta = (self.Hc*self.Rc**2/2)* (-((1/a-1)*a**2*cos(theta))/(sqrt(1-(1/a-1)**2*(sin(theta))**2)) - a**2 + (1-a)**2*(-cos(2*theta))-a*(1-a)*cos(theta)*sqrt(1-(1/a-1)**2*(sin(theta))**2)+ ((1/a-1)**2*a*(1-a)*(sin(theta))**2*cos(theta)  )/(sqrt(1-(1/a-1)**2*(sin(theta))**2)) + 1   ) - (self.b*self.Hc*self.Rc/2)*( (1-a)*sin(theta)+ ((1-a)**2*sin(theta)*cos(theta))/(sqrt((1-a)**2*(cos(theta))**2+2*a-1))  )
        
        if V < 0.:
            raise Exception('V negative %f' % (V))
            V = 0.0
        return V,dVdtheta
        
    def Vs_dVs_ParkYC(self,theta):
        #Suction chamber
        if theta >= 0.0 and theta <= 2*pi:    
            V,dVdtheta = self.V_dV_ParkYC(theta)
            return abs(V),abs(dVdtheta)
        elif theta > 2*pi and theta <= 4*pi:
            V,dVdtheta = self.V_dV_ParkYC(4*pi - theta)
            return abs(V),-abs(dVdtheta)
        
    def Vc_dVc_ParkYC(self,theta):
        #Compression chamber
        if theta >= 0.0 and theta <= 2*pi:    
            V,dVdtheta = self.V_dV_ParkYC(2*pi - theta)
            return abs(V),-abs(dVdtheta)
        elif theta > 2*pi and theta <= 4*pi:
            V,dVdtheta = self.V_dV_ParkYC(theta-2*pi)
            return abs(V),abs(dVdtheta)

    def Vs_dVs_ParkYC_2stage(self,theta):
        #Second cylinder suction chamber
        if theta >= 0.0 and theta <= pi:    
            V,dVdtheta = self.V_dV_ParkYC(theta + pi)
            return abs(V),abs(dVdtheta)        
        elif theta > pi:
            V,dVdtheta = self.V_dV_ParkYC(theta - pi)
            return abs(V),abs(dVdtheta) 

    def Vc_dVc_ParkYC_2stage(self,theta):
        #Second cylinder compression chamber
        if theta >= 0.0 and theta <= pi:    
            V,dVdtheta = self.V_dV_ParkYC(2*pi - theta - pi)
            return abs(V),-abs(dVdtheta)        
        elif theta > pi:
            V,dVdtheta = self.V_dV_ParkYC(2*pi - theta + pi)
            return abs(V),-abs(dVdtheta)         
        
    def x_theta(self,x):
        """
        Vane displacement with the crankshaft angle
        """
        a = self.Rr/self.Rc
        return self.Rc*(1-(1-a)*cos(x) - sqrt((1-a)**2*(cos(x))**2+2*a-1 ) )
 
    def xv_dot(self,theta, stage=1):
        """
        Function for calculating velocity of the vane.
        """
        R_r = self.Rr
        R_v = self.rv
        e = self.e
        z = 1 - (e*sin(theta)/(R_r+R_v))**2
        alphadot = e*cos(theta)*self.omega/(sqrt(z)*(R_r+R_v))
        
        return e*sin(theta)*(alphadot+self.omega)

    def V_disp(self,alpha_s,alpha_d):
        """
        Zheng et al. 'Experimental verification of a rolling-piston expander that
        applied for low-temperature ORC', Applied Energy 112(2013)1265-1274.
        
        ParkY.C. 'Transient analysis of a variable speed rotary compressor' 
        Energy Conversion and Management 51(2010) 277-287 
        
        """
        V_theta_s = self.Vs_dVs_ParkYC(alpha_s)[0]
        V_theta_d = self.Vc_dVc_ParkYC(2*pi - alpha_d)[0]
        Int_x_dtheta, err = quad(self.x_theta,0,2*pi)
        Vdisp = pi*self.Hc*(self.Rc**2 - self.Rr**2) - self.b*self.Hc*Int_x_dtheta - (V_theta_s + V_theta_d)

        if hasattr(self.mech,'cylinder_config') and self.mech.cylinder_config == 'dual-cylinder':
            return 2*Vdisp
        elif hasattr(self.mech,'cylinder_config') and self.mech.cylinder_config == 'single-cylinder':
            return Vdisp
        else:
            return Vdisp

    def V_shell(self, theta):
        #Approximation of shell volume
        V_shell_tot = pi*self.Rshell**2*self.Hshell
        Vcylinder = pi*self.Rc_outer**2*self.Hc
        Vmotor = pi*self.Rshell**2*self.Hc/3
        V = (V_shell_tot - Vcylinder - Vmotor)
        dV = 0.0
        return V,dV

    def xv_simplified(self,theta, stage = 1):
        """
        Distance that the vane extends into the cylinder
        Assumption: rv = 0.0
        """
        if stage == 2:
            R_c = self.Rc_2
            R_r = self.Rr_2
            e = self.e_2
        elif stage == 1:
            R_c = self.Rc
            R_r = self.Rr
            e = self.e
        
        return R_c - R_r * cos(asin(e / R_r * sin(theta))) - E * cos(theta)
    
    def A_vb(self,theta, geo):
        """
        Leakage area through vertical clearance between vane and cylinder
        """
        return 2*self.delta_vb*self.x_theta(theta,geo)

    def A_vt(self):
        """
        Leakage area through radial clearance between vane tip and roller
        """
        return self.hv*self.delta_vt

    def A_f(self):
        """
        Leakage area between compression chamber and muffler
        """
        return pi/4*self.D_muff**2
        
    def A_lc(self):
        """
        Leakage area between upper and lower shell
        """
        return self.A_ls

    def A_uc(self):
        """
        Leakage area from discharge port
        """
        return pi/4*self.D_uc**2

    def A_rc(self,theta, delta):
        """
        Function for calculating area of leakage path across roller.
        """
        x = self.e*sin(theta)/self.Rr
        alpha = asin(x)
        beta = alpha + theta
        
        return 2*(2*pi - beta)*self.Rr_i*delta
    
    def A_10(self):
        """
        Leakage area between accumulator and suction pipe
        """
        return self.A_sp

    def A_21(self,theta, stage=1):
        """
        Leakage area between suction chamber and suction pipe
        """
        l_theta_max = pi*self.D_sp
        
        if theta <= self.theta_s1:
            l_theta = 0.0
            A21 = 0.0
        elif theta > self.theta_s1 and theta <= self.theta_sm:
            l_theta = l_theta_max/2*(theta - self.theta_s1)/(self.theta_sm - self.theta_s1)
            phi = theta - self.theta_s1
            z_theta = self.xv_simplified(phi, stage)
            A21 = z_theta*l_theta/2
        elif theta > self.theta_sm and theta <= self.theta_s2:
            l_theta = l_theta_max*(1 - 0.5*(self.theta_s2 - theta)/(self.theta_sm - self.theta_s1))
            phi = theta - self.theta_s1
            z_theta = self.xv_simplified(phi, stage)
            A21 = z_theta*l_theta/2
        elif theta > self.theta_s2:
            l_theta = l_theta_max
            phi1 = theta - self.theta_s1
            phi2 = theta - self.theta_s2
            z_theta = (self.xv_simplified(phi1, stage) + self.xv_simplified(phi2, stage))/2
            A21 = z_theta*l_theta

        if A21 > self.A_sp: 
            A21 = self.A_sp
        
        return A21

    def A_31(self,theta, stage=1):
        """
        Leakage area between compression chamber and suction pipe
        """
        l_theta_max = pi*self.D_sp
    
        if theta <= self.theta_s1:
            l_theta = l_theta_max
        elif theta > self.theta_s1 and theta <= self.theta_sm:
            l_theta = l_theta_max*(1 - 0.5*(theta-self.theta_s1)/(self.theta_sm-self.theta_s1))
        elif theta > self.theta_sm and theta < self.theta_s2:
            l_theta = l_theta_max/2*(self.theta_s2 - theta)/(self.theta_sm - self.theta_s1)
        elif theta >= self.theta_s2:
            l_theta = 0
    
        phi = self.theta_s2 - theta
        z_theta = self.xv_simplified(phi, stage)
        A31 = z_theta*l_theta/2
        
        if A31 > self.A_sp: 
            A31 = self.A_sp
        
        return A31

    def A_32(self,theta, stage=1):
        """
        Flank leakage area between compression chamber and suction chamber
        """
        if theta <= self.theta_d1:
            H32 = self.Hc
        elif theta > self.theta_d1 and theta < self.theta_dm:
            H32 = self.Hc - self.hp*(theta-self.theta_d1)/(self.theta_dm-self.theta_d1)
        else:
            H32 = self.Hc - self.hp*(self.theta_d2-theta)/(self.theta_dm-self.theta_d1)
        
        if H32 > self.Hc: 
            H32 = self.Hc
        
        return H32*self.delta_gap

    def A_42_ref(self,stage=1):
        """
        Reference value used to calculate flow coefficient from clearance volume to suction chamber
        """
    
        phi = self.theta_d2 - self.theta_d1
        z_theta = self.xv_simplified(phi, stage)

        return self.l_dm*z_theta/2

    def A_42(self,theta, stage=1):
        """
        Leakage area between clearance volume and suction chamber
        """
        phi = self.theta_d2 - self.theta_d1
        z_theta = self.xv_simplified(theta, stage)
    
        if theta < self.theta_d1:
            A42_side = 0.0
            A42_t4 = 0.0
            A42 = 0.0
        elif theta >= self.theta_d1 and theta < self.theta_d2:
            l_theta = self.l_dm * (self.theta - self.theta_d1)/(self.theta_d2 - self.theta_d1)
            phi = theta - self.theta_d1
            z_theta = self.xv_simplified(phi, stage)
            A42_side = l_theta*z_theta/2
            A42_t4 = z_theta*(self.D_dp- self.D_dpc)/2
            A42 = A42_side + A42_t4
        elif theta >= self.theta_d2:
            phi1 = self.theta_d2 - self.theta_d1
            phi2 = theta - self.theta_d1
            z_theta1 = self.xv_simplified(phi1, stage)
            z_theta2 = self.xv_simplified(phi2, stage)
            A42_t4 = z_theta2*(self.D_dp - self.D_dpc)/2
            A_42 = self.l_dm*z_theta1/2 # + A42_t4

        return A42    

    def A_43_ref(self,theta, stage=1):
        """
        Reference value used to calculate flow coefficient from clearance volume to compression chamber
        """
        A43a_max = pi/8*self.D_dp**2
        phi1 = self.theta_d1 - theta
        phi2 = self.theta_d2 - theta
        phim = self.theta_dm - theta
        
        z_thetam = self.xv_simplified(phim, stage)
        if z_thetam > self.D_dp/2:
            z_thetam = self.D_dp/2
        
        if theta < self.theta_d1:
            z_theta2 = (self.xv_simplified(phi1, stage) + self.xv_simplified(phi2, stage))/2
        else:
            z_theta2 = self.xv_simplified(phi2, stage)/2
        
        A43_side = z_thetam*self.hp/4
        A43_a = self.D_dp*z_theta2
        if A43_a > A43a_max:
            A43_a = A43a_max
        A43_b = pi*self.D_dp*self.hp/4
        A43_ref = sqrt(A43_a**2 + A43_b**2) + A43_side
        
        return A43_ref
       
    def A_43(self,theta, stage=1):
        """
        Leakage area between clearance volume and compression chamber
        """
        A43a_max = pi/8*self.D_dp**2
        phi1 = self.theta_d1 - theta
        phi2 = self.theta_d2 - theta
        phim = self.theta_dm - theta
        
        z_thetam = self.xv_simplified(phim, stage)
        if z_thetam > self.D_dp/2:
            z_thetam = self.D_dp/2
        
        if theta < self.theta_d1:
            z_theta2 = (self.xv_simplified(phi1, stage) + self.xv_simplified(phi2, stage))/2
        else:
            z_theta2 = self.xv_simplified(phi2, stage)/2

        l_theta = self.l_dm*(self.theta_d2 - theta)/(self.theta_d2 - self.theta_d1)
        if l_theta < 0.0:
            l_theta = 0.0
        elif l_theta > self.l_dm:
            l_theta = self.l_dm
        
        A43_side = z_thetam*self.hp/4
        A43_a = self.D_dp*z_theta2
        if A43_a > A43a_max:
             A43_a = A43a_max
        
        A43_b = pi*self.D_dp*self.hp/4
        A43_c1 = sqrt(A43_a**2 + A43_b**2) + A43_side
        A43_c2 = l_theta*z_theta2
        
        if A43_c1 > A43_c2:
            A43 = A43_c2
        else:
            A43 = A43_c1
        
        return A43
    
    def calculate_heat_transfer_area(self,theta):
        """
        Calculation of chamber total heat transfer area
        """
        a = self.Rr/self.Rc
        hv_theta = self.Rc*(1-(1-a)*cos(theta)-sqrt((1-a)**2*(cos(theta))**2) +2*a-1)
        f_theta = (1-a**2)*theta-0.5*(1-a)**2*sin(2*theta) - a**2*asin((1/a-1)*sin(theta)) - a*(1-a)*sin(theta)*sqrt(1-(1/a-1)**2*(sin(theta))**2)
        phi = pi/2 + theta - asin((self.Rc - hv_theta + self.e*cos(pi - theta))/self.Rr )
        
        #delta = self.Rc- (2*self.e*cos(theta-phi)+sqrt(4*e**2*(cos(theta-phi))**2-4*(e**2 - self.Rr**2)))/2
        
        A = theta*(self.Rc + self.Rr)*self.Hc+ 2*0.5*self.Rc**2*f_theta
        
        if A < 0.0:
            raise Exception('A_ht < 0.0')
        
        return A

    def A_HT_s(self,theta):
        #Suction chamber
        if theta >= 0.0 and theta <= 2*pi:    
            A = self.calculate_heat_transfer_area(theta)
            return A
        elif theta > 2*pi and theta <= 4*pi:
            A = self.calculate_heat_transfer_area(4*pi-theta)
            return A
 
    def A_HT_c(self,theta):
        #Compression chamber
        if theta >= 0.0 and theta <= 2*pi:    
            A = self.calculate_heat_transfer_area(2*pi - theta)
            return A
        elif theta > 2*pi and theta <= 4*pi:
            A = self.calculate_heat_transfer_area(theta-2*pi)
            return A  

    def A_HT_s_2stage(self,theta):
        #Second cylinder suction chamber
        if theta >= 0.0 and theta <= pi:    
            A = self.calculate_heat_transfer_area(theta + pi)
            return A
        elif theta > pi:
            A = self.calculate_heat_transfer_area(theta - pi)
            return A
 
    def A_HT_c_2stage(self,theta):
        #Second cylinder compression chamber
        if theta >= 0.0 and theta <= pi:    
            A = self.calculate_heat_transfer_area(2*pi - theta - pi)
            return A
        elif theta > pi:
            A = self.calculate_heat_transfer_area(2*pi - theta + pi)
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

    def alp(self,p,T):   
        """
        Solubility Refrigerant in oil
        p: pressure [kPa]
        T: temperature [K]
        """                                  
        if self.Ref == 'R410a' or self.Ref == 'HEOS::R410a' or self.Ref == 'REFPROP::R401a.mix':
            #Solubility of the R410A and PVE Oil
            P0 = p*1e-3  #pressure [MPa]
            T0 = T - 273.15 #Temperaure [C]
            
            Sol20 = 3.49058 + 14.33071*P0 + 3.94294*P0**2
            Sol40 = 5.57178 + 0.76551*P0 + 5.97508*P0**2
            Sol60 = 4.81862 - 0.34116*P0 + 3.51925*P0**2
            Sol80 = 2.88481 + 2.19068*P0 + 1.54295*P0**2
            Sol100 = 2.2739 + 2.24531*P0 + 1.07309*P0**2
            Sol120 = 0.56494 + 3.15252*P0 + 0.65613*P0**2
            
            if T0 < 20:
                alp = Sol20
            elif 20 <= T0 and T0 < 40:
                alp = Sol40 - (Sol40 - Sol20)*(40 - T0)/(40 - 20)
            elif T0 >= 40 and T0 < 60:
                alp = Sol60 - (Sol60 - Sol40)*(60 - T0)/(60 - 40)
            elif T0 >= 60 and T0 < 80:
                alp = Sol80 - (Sol80 - Sol60)*(80 - T0)/(80 - 60)
            elif T0 >= 80 and T0 < 100:
                alp = Sol100 - (Sol100 - Sol80)*(100 - T0)/(100 - 80)
            elif T0 >= 100 and T0 <= 120:
                alp = Sol120 - (Sol120 - Sol100)*(120 - T0)/(120 - 100)
            elif T0 > 120:
                alp = Sol120

            if alp <= 0:
                alp = 0
            if alp > 10:
                alp = 10
    
        else:
            alp = 0.0
        
        return alp/100

    def Oil_viscosity(self,p,T):
        """
        Used to calculate oil viscosity as a function of temperature and pressure.
        p : pressure kPa
        T : temperature K
        """
        alp1 = self.alp(p, T)*100                                                                            
        T1 = T - 273.15                                                                                     
        
        if T1 < 60:                                                                                                                         
            #0<=T< 60
            mu_10 = 196.05-13.05438*T1 + 0.39545*T1**2 - 0.00617*T1**3 + 0.0000478516*T1**4 - 0.000000145573*T1**5
        elif T1 >= 60:
            #60<=T<=100
            mu_10 = 59.9 - 1.52308*T1 + 0.01471*T1**2 - 0.0000504167*T1**3

        if T1 < 20:                                                                                      
            #-20<=T<20
            mu_20 = 44.75 - 2.83417*T1 + 0.10098*T1**2 - 0.00198*T1**3 + 0.0000212222*T1**4 - 0.000000116146*T1**5 + 0.000000000256076*T1**6
        elif T1 >= 20:                                                                                  
            #20<=T<=100
            mu_20 = 24.682 - 0.5437*T1 + 0.00472*T1**2 - 0.0000145833*T1**3
        
        if T1 < 60:                                                                                      
            #-20<=T<60
            mu_30 = 13.78 - 0.69392*T1 + 0.03069*T1**2 - 0.000867812*T1**3 + 0.0000135104*T1**4 - 0.000000105937*T1**5 + 0.000000000326823*T1**6
        elif T1 >= 60:                                                                                 
            #60<=T<=100
            mu_30 = 10.15 - 0.184*T1 + 0.00128*T1**2 - 0.0000025*T1**3

        mu_40 = 6.60955 - 0.20951*T1 + 0.00385*T1**2 - 0.0000303243*T1**3 + 0.0000000262784*T1**4 + 0.000000000507812*T1**5
        mu_50 = 3.18206 - 0.07856*T1 + 0.00203*T1**2 - 0.0000350712*T1**3 + 0.000000323208*T1**4 - 0.00000000115685*T1**5
        
        if alp1 < 20:
            mu_oil = (mu_20 - ((20 - alp1)*(mu_20 - mu_10)/10))*0.001
        elif alp1 >= 20 and alp1 < 30:
            mu_oil = (mu_30 - ((30 - alp1)*(mu_30 - mu_20)/10))*0.001
        elif alp1 >= 30 and alp1 < 40:
            mu_oil = (mu_40 - ((40 - alp1)*(mu_40 - mu_30)/10))*0.001
        elif alp1 >= 40:                                                                                   
            #And alp1 <= 90 Then
            mu_oil = (mu_50 - ((50 - alp1)*(mu_50 - mu_40)/10))*0.001
        else:
            print(alp1, T1)
            raise Exception('Out of range')
        
        return mu_oil

    def Inlet(self, FlowPath, **kwargs):
        #if self.theta >= self.psi:
        FlowPath.A=self.A_suction
        mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
        return mdot

    def Outlet(self, FlowPath, **kwargs):
        FlowPath.A=self.A_discharge
        mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
        return mdot

    def Suction(self,FlowPath,**kwargs):
        if FlowPath.key_up=='A':
            # pressure in compressor higher than the inlet line
            # valve is closed - no flow
            return 0.0
        if self.theta >= self.theta_s1:
            FlowPath.A=self.A_suction
            mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
            return mdot
        else:
            return 0.0

    def Suction_2stage(self,FlowPath,**kwargs):
        if FlowPath.key_up=='C':
            # pressure in compressor higher than the inlet line
            # valve is closed - no flow
            return 0.0
        if self.theta <= pi or self.theta >= self.theta_s1+pi:
            FlowPath.A=self.A_suction
            mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
            return mdot
        else:
            return 0.0
        
    def Discharge(self,FlowPath,**kwargs):
        if FlowPath.key_down=='B':
            # pressure in compressor lower than the discharge line
            # valve is closed - no flow
            return 0.0
        else:
            if self.theta <= self.theta_d2:
                FlowPath.A=self.discharge_valve.A()
                mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
                return mdot
            else:
                 return 0.0

    def Discharge_2stage(self,FlowPath,**kwargs):
        if FlowPath.key_down=='D':
            # pressure in compressor lower than the discharge line
            # valve is closed - no flow
            return 0.0
        else:
            if self.theta <= self.theta_d2 - pi:    
                FlowPath.A=self.discharge_valve_upper.A()
                mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)                
                return mdot
            else:
                return 0.0

    def CompressibleGasLeakage_vb(self, FlowPath, **kwargs):
        """
        From compression to suction chamber through vertical clearance between vane and cylinder
        """
        Cv = self.Cflow_vb
        FlowPath.A=self.delta_vb*self.b
        mdot=Cv*flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
        return mdot

    def CompressibleGasLeakage_vt(self, FlowPath, **kwargs):
        """
        From compression to suction chamber through radial clearance between vane tip and rolling piston
        """        
        Cv = self.Cflow_vt
        FlowPath.A=self.delta_vt*self.Hc
        mdot=Cv*flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
        return mdot

    def CompressibleGasLeakage_32(self, FlowPath, **kwargs):
        """
        From compression to suction chamber through radial clearance between cylinder wall and 
        rolling piston
        """        
        Cv = self.Cflow_32
        FlowPath.A=self.A_32(self.theta)
        mdot=Cv*flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)        
        return mdot

    def CouettePoiseuilleLeakage_VaneSlot(self, FlowPath, **kwargs):        
        State_up = FlowPath.State_up
        State_down = FlowPath.State_down        
        
        T_up = State_up.get_T()
        P_up = State_up.get_p()
        rho_up = State_up.get_rho()
        
        T_down = State_down.get_T()
        P_down = State_down.get_p()  
        rho_down = State_up.get_rho()      
        
        if FlowPath.key_up=='shell':
            # Vane slot chamber
            x_up = self.x_slot            
            if FlowPath.key_down=='A':
                # Suction chamber
                x_down = self.xs 
            elif FlowPath.key_down=='B':
                # Compression chamber
                x_down = self.xc        
        elif FlowPath.key_up=='B':
            # Compression chamber
            x_up = self.xc
            # Vane slot chamber
            x_down = self.x_slot
        elif FlowPath.key_up=='A':
            # Suction chamber
            x_up = self.xs
            # Vane slot chamber
            x_down = self.x_slot        
        else:
            raise Exception('VaneSlotLeakage condition not defined')
        
        #Volume flow rate calculated using mixed Plane Couette and Poiseuille flow.
        mu_gas_up = State_up.get_visc()
        mu_oil_up = self.Oil_viscosity(P_up, T_up)
        mu_up = x_up*mu_oil_up + (1 - x_up)*mu_gas_up
        mu_gas_down = State_down.get_visc()
        mu_oil_down = self.Oil_viscosity(P_down, T_down)
        mu_down = x_down * mu_oil_down + (1 - x_down) * mu_gas_down
        mu_mean = (mu_up + mu_down)/2
        
        xv_dot = self.xv_dot(self.theta)
        
        Vdot_eb = self.hv*self.delta_vb*(xv_dot/2 + (P_up - P_down)*self.delta_vb**2/(12*mu_mean*self.L_slot))
        
        if Vdot_eb >= 0:
            rho_mix = 1/(x_up/self.rho_oil + (1 - x_up)/rho_up)
        else:
            rho_mix = 1/(x_down/self.rho_oil + (1 - x_down)/rho_down)
        
        mdot = Vdot_eb*rho_mix

        return mdot
        
    def LaminarViscousLeakage_Roller(self, FlowPath, **kwargs):
        """
        Used to calculate mass flow rate for leakage across the roller, modeled as viscous flow.
        """
        # Rotation angle
        theta = self.theta
        
        State_up = FlowPath.State_up
        State_down = FlowPath.State_down
        
        T_up = State_up.get_T() #[K]
        P_up = State_up.get_p() #[kPa]
        rho_up = State_up.get_rho() #[kg/m3]

        T_down = State_down.get_T() #[K]
        P_down = State_down.get_p() #[kPa] 
        rho_down = State_up.get_rho() #[kg/m3]

        if FlowPath.key_up=='shell':
            # Oil concentration across roller
            x_up = self.x_roller
            if FlowPath.key_down=='A':
                # Suction chamber
                x_down = self.xs 
            elif FlowPath.key_down=='B':
                # Compression chamber
                x_down = self.xc        
        elif FlowPath.key_up=='B':
            # Compression chamber
            x_up = self.xc
            # Oil concentration across roller
            x_down = self.x_roller
        elif FlowPath.key_up=='A':
            # Suction chamber
            x_up = self.xs
            # Oil concentration across roller
            x_down = self.x_roller        
        else:
            raise Exception('RollerLeakage condition not defined')        
        
        mu_gas_up = State_up.get_visc() #[Pa-s]
        mu_oil_up = self.Oil_viscosity(P_up, T_up)
        mu_up = x_up*mu_oil_up + (1 - x_up)*mu_gas_up
        mu_gas_down = State_down.get_visc()
        mu_oil_down = self.Oil_viscosity(P_down, T_down)
        mu_down = x_down * mu_oil_down + (1 - x_down) * mu_gas_down
        mu_mean = (mu_up + mu_down)/2
        
        x = self.e*sin(theta)/self.Rr
        alpha = asin(x)
        beta = alpha + theta
        
        if FlowPath.key_down=='B' or FlowPath.key_up=='B' :
            A = 2*(2*pi - beta)*(self.Rr + self.Rr_i)/2*self.delta_side
        elif FlowPath.key_up=='A' or FlowPath.key_down=='A':
            A = 2*beta*(self.Rr + self.Rr_i)/2*self.delta_side
        
        if FlowPath.key_up=='shell' and FlowPath.key_down=='B' :
            DeltaP = (P_up - P_down) # [kPa]
            Vdot_rc = 2*pi*self.delta_side**3*DeltaP/(6*mu_mean*log(self.Rr/self.Rr_i))*((2*pi - beta)/(2 * pi))
            rho_mix = 1/(x_up/self.rho_oil + (1 - x_up)/rho_up)
            mdot = rho_mix*Vdot_rc
        else:#FlowPath.key_up=='B' or FlowPath.key_up=='A':            
            Cv = self.Cflow_roller
            mdot=Cv*flow_models.IsentropicNozzle(A,FlowPath.State_up,FlowPath.State_down)
        # else:
        #     raise Exception('RollerLeakage condition not defined')
        return mdot        
              
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
                    V,dVdtheta = self.Vs_dVs_ParkYC(theta)
                    A_ht = self.A_HT_s(theta) #[m2]
                    
                elif key == 'B':
                    #Compression chamber
                    V,dVdtheta = self.Vc_dVc_ParkYC(theta)
                    A_ht = self.A_HT_c(theta) #[m2]

                if key == 'C':
                    #Suction chamber
                    V,dVdtheta = self.Vs_dVs_ParkYC_2stage(theta)
                    A_ht = self.A_HT_s_2stage(theta) #[m2]
                    
                elif key == 'D':
                    #Compression chamber
                    V,dVdtheta = self.Vc_dVc_ParkYC_2stage(theta)
                    A_ht = self.A_HT_c_2stage(theta) #[m2]
                
                #Assign lumped temperature to wall
                T_w  = self.Tlumps[0] #[K]
                
                #Calculate thermophysical properties
                cp = self.CVs[key].State.cp # kJ/kg/K    
                rho = self.CVs[key].State.rho #[kg/m3]
                k = self.CVs[key].State.k #[kW/m-K]
                mu = self.CVs[key].State.visc #[Pa-s]
                T = self.CVs[key].State.T #[K]
                
                #Avoid division by zero when Aht or Dh are zero
                r_avg =  (self.Rc + self.Rr)/2  # self.r_chamber_mean(theta)              
                if A_ht > 0.0 and V > 0.0:
                    D_h = 4*V/A_ht # Liu (1994) Eq. 3.4.2 [m]                    
                    if D_h < 0.0:
                        raise Exception('D_h < 0.0')
                    
                    if D_h > 0.0:
                        u = D_h*self.omega/2  # Liu (1994) Eq. 3.4.3 #[m/s]
                        Pr = (cp*mu)/k #[-]
                        Re = (rho*u*D_h)/mu #[kg/m3*m/s*m/Pa/s]=[-]
    
                        h_c = 0.025*(k/D_h)*Pr**(0.4)*Re**(0.8)*(1.0+1.77*(D_h/r_avg)) #[kW/m2/K] Liu (1994) Eq. 3.4.8
                        Q = h_c*A_ht*(T_w-T)   #Net heat into control volume [kW]
                    else:
                        Q = 0.0 #kW
                else:
                    Q = 0.0 #kW

                if key == 'shell':                    
                    T_w = self.Tlumps[0] # [K]
                    T = self.CVs[key].State.T # [K]
                    P_film = self.CVs[key].State.p # [kPa]
                    
                    # Vertical side shell
                    h_c_vertical = HTC('vertical_plate',T_w,T,P_film,self.Ref,self.Hshell)
                    A_shell_vertical = 2*pi*self.Rshell*self.Hshell # m2
                    Q_vertical = h_c_vertical*A_shell_vertical*(T_w-T) 

                    # Top Surface
                    h_c_topshell = HTC('horizontal_plate',T_w,T,P_film,self.Ref,self.Rshell*2,PlateNum='upper_surface')
                    A_shell_top = pi*self.Rshell**2 # m2
                    Q_topshell = h_c_topshell*A_shell_top*(T_w-T)                     

                    # Bottom Surface
                    h_c_bottomshell = HTC('horizontal_plate',T_w,T,P_film,self.Ref,self.Rshell*2,PlateNum='lower_surface')
                    A_shell_bottom = pi*self.Rshell**2 # m2
                    Q_bottomshell = h_c_bottomshell*A_shell_bottom*(T_w-T)  
                    
                    # Total heat transfer
                    Q = Q_vertical + Q_topshell + Q_bottomshell
                                    
                Qarray.append(Q)

            return arraym(Qarray)
        
        elif self.HT_on == False:
            return arraym([0.0]*self.CVs.N)
            
        else:
            raise Exception('HT_on is either True or False')

        
    def calculate_force_terms(self,p_shell = None):
        """
        Calculate the force profiles, mean forces, moments, etc.
        
        Parameters
        ----------
        p0 : float, or class instance
            If a class instance, must provide a function __call__ that takes as its first input the Rolling class
        
        """
        if not isinstance(p_shell,float):
            raise ValueError('calculate_force_terms must get a float top vane pressure for now')
        
        Nstep = self.Itheta+1
        _slice = range(self.Itheta+1)
        theta = self.t[_slice]    
        
        e = self.e
        omega = self.omega
        epsilon = e/self.Rr
        hv = self.hv         #vane height 
        b = self.b          #vane width
        rv = self.rv
        Rr = self.Rr
        Re = self.Re
        Rs = self.Rs
        Hc = self.Hc
        rho_v = 7200 #kg/m3
        rho_r = 7200 #kg/m3
        I_r = pi*rho_r*(Rr**4 - Re**4)*Hc/2 #Moment of inertia of rolling piston
        
        #Volume of the vane from Liu PhD Thesis
        av = rv - sqrt(rv**2 - b**2/4 )
        theta_v = asin(b/(2*rv))
        A5 = rv**2*theta_v- b/2*sqrt(rv**2 - b**2/4)
        Av = b*(hv-av)
        Vv = Hc*(Av + A5)
        mv = rho_v*Vv

        #vane displacement
        xv = e*(1-np.cos(theta)) + Rr*(1 - np.sqrt(1-epsilon**2*(np.sin(theta))**2))
        #vane sliding velocity
        xdot_v = e*omega*(np.sin(theta) + (epsilon*np.sin(2*theta))/(2*np.sqrt(1-epsilon**2*(np.sin(theta))**2))) 
        #vane acceleration
        self.xdotdot_v = e*omega**2*(np.cos(theta) + (epsilon*np.cos(2*theta))/(np.sqrt(1- epsilon**2*(np.sin(theta))**2)) + (epsilon**3*np.sin(2*theta))/(4*(np.sqrt( (1-epsilon**2*(np.sin(theta))**2)**3  ))) )
        alpha = np.arcsin(e*np.sin(theta)/(Rr+rv))
        pc = self.p[1][_slice]*1000 #Pa
        pe = self.p[0][_slice]*1000 #Pa

        Tc = self.T[1][_slice] #K
        Te = self.T[0][_slice] #K

        #Forces on the vane
        self.Fc = Hc*(b*p_shell - pe*(b/2 - rv*np.sin(alpha)) - pc*(b/2 + rv*np.sin(alpha)))
        self.Fh = Hc*xv*(pc-pe)
        self.F_lv = -mv*self.xdotdot_v
        x0 = 0.0 #initial extension spring
        K_v = 1150    #spring stiffness [N/m] can be obtained from Ooi 1994 800-3000 N/m
        self.Fk = K_v*(x0 - np.fabs(xv))
        
        #Forces on the roller
        m_r = rho_r*(pi*(Rr**2-Re**2)*Hc)
        #Centrifugal force on roller
        self.F_lp = m_r*e*omega**2

        """
        From Yanagisawa and Shimizu "Friction losses in rolling piston type compressor III"        
        mu_v = 0.15 - 35*np.sqrt(mu_v/(Fn/Hc))
        """
        mu_v = 0.16  
        mu_s = 0.18 #typically 0.16-0.22
        delta_rv = rv*(1-np.cos(alpha))
        self.Fn = (mu_s*self.Fh*(hv+mu_s*b)+(self.Fk+self.Fc+self.F_lv)*(xv-hv))/( mu_s*(np.sin(alpha)-mu_v*np.cos(alpha))*(xv+hv+mu_s*b-2*delta_rv) + (np.cos(alpha) + mu_v*np.sin(alpha))*(xv - hv - 2*mu_s*rv*np.sin(alpha))        )

        self.F_R1 = 0.5*self.Fh + 1/(2*mu_s)*(self.Fk + self.Fc + self.F_lv) + 0.5*self.Fn*(mu_v*np.cos(alpha) - (mu_v/mu_s)*np.sin(alpha) - np.sin(alpha) - np.cos(alpha)/mu_s)
        self.F_R2 = -0.5*self.Fh + 1/(2*mu_s)*(self.Fk + self.Fc + self.F_lv) - 0.5*self.Fn*(mu_v*np.cos(alpha) + (mu_v/mu_s)*np.sin(alpha) - np.sin(alpha) + np.cos(alpha)/mu_s)
        
        #Friction force generated by radial force of vane on roller
        self.Ft = mu_v*self.Fn
        self.Ft_R1 = mu_s*self.F_R1
        self.Ft_R2 = mu_s*self.F_R2
        # force caused by pressure difference between suction and discharge chambers
        self.Fg = 2*Rr*Hc*(pc-pe)*np.sin((theta+alpha)/2)
        #Angular position of Fg
        self.theta_Fg = (theta - alpha)/2
        
        #Force balance equations on roller
        self.Fr = self.Fg*np.cos((theta+alpha)/2) - self.Fn*np.cos(theta+alpha) - self.Ft*np.sin(theta+alpha) + self.F_lp
        self.Ftheta = -self.Fg*np.sin((theta+alpha)/2) + self.Fn*np.sin(theta+alpha) - self.Ft*np.cos(theta+alpha)
        
        #Resultant
        self.F = np.sqrt(self.Fr**2 + self.Ftheta**2)
        self.theta_f = theta + np.arctan(self.Ftheta/self.Fr)
        
        #Gas torque
        self.Mg = e*self.Fg*np.sin((theta+alpha)/2)
        #Torque generated by radial force on the center of roller imposed by vane 
        self.Mn = -e*self.Fn*np.sin(theta+alpha)
        #Friction torque generated by mu_v*Fn = Ft
        #TODO: check sign
        self.Mn_t = e*self.Ft*np.cos(theta+alpha)
        #Moment acting on the crankshaft by the roller
        self.Mf = self.Mg + self.Mn + self.Mn_t

        #Assming steady-state conditions
        Toil = (Tc+Te)/2 #K
        rho_oil = -0.0034*(Toil-273.15)**2-0.0776*(Toil-273.15) + 960.88
        #Vogel equation
        mu_oil = rho_oil*(1.742e-4)*np.exp(3664/(Toil-10.27))*1e-6  #Toil [K]
        #Clearences
        delta_e = 0.25e-3 #[m]
        delta_1 = 0.25e-3 #[m]
        delta_2 = 0.25e-3 #[m]
        delta_3 = 0.25e-3 #[m]
        delta_4 = 0.25e-3 #[m]

        #Simplified analysis roller rotational speed
        l_r = np.sqrt(Re**2+e**2-2*Re*e*np.cos(pi-self.theta_f+theta))
        cos_thetar = np.sqrt(1 - (e/l_r*np.sin(self.theta_f - theta))**2)
        self.omega_r = (2*pi*mu_oil*omega*Hc*Re**3/delta_2 - Rr*self.Ft)*delta_2*delta_e/(2*pi*mu_oil*(Hc*Re**3 + (Hc*Re**3*delta_e)*delta_2))
        
#         self.omega_r = (mu_oil*omega*Hc*Re**3/delta_2 - Rr*self.Ft/(2*pi))/(mu_oil*Hc*Re**3/delta_2 + mu_oil*(Rr**4 - Re**4)/delta_1)

        #Relative angular velocity
        DELTA_omega = omega - self.omega_r
        #Frictional loss between eccentric and the inner surface of the roller Ler
        self.Ler = DELTA_omega*(2*pi*mu_oil*DELTA_omega)*Hc*Re**3/delta_2
        
        #Frictional loss between the roller face and the cylinder head faces Lrc
        self.Lrc = omega*(2*pi*self.omega_r*mu_oil*(Rr**4-Re**4))/delta_1
        
        #Frictional loss between eccentric face and the cylinder head faces Lec
        self.Lec = omega*(pi*omega*mu_oil*(2*Re**4-Rs**4))/delta_3
        
        #Frictional loss between vane tip and roller Lv
        self.Lv = self.Ft*(Rr*self.omega_r + (e*omega*np.cos(theta))/np.cos(alpha))
        
        #Frictional loss between vane sides and vane slot Ls
        alphadot = 1/omega*(e*np.cos(theta)/(Rr+rv))/(np.sqrt(1 - (e*np.sin(theta)/(Rr+rv))**2))
        xdot = e*omega*np.sin(theta) + (Rr+rv)*alphadot*np.sin(alpha)
        self.Ls = (self.Ft_R1 + self.Ft_R2)*np.fabs(xdot)
        #Frictional loss between shaft and upper journal bearing  Lb2
        h_mj2 = 0.067 # length bearing[m]
        self.Lbs = omega*(2*pi*mu_oil*omega*h_mj2*Rs**3)/delta_4
        
        #Frictional loss between shaft and lower bearing , Lb1
        h_sj = 0.025 # length bearing[m]
        self.Lbl = omega*(2*pi*mu_oil*omega*h_sj*Rs**3)/delta_4

        #Lower and upper balancer distance and eccentricity
        l_b1 = 0.1241 #[m]
        l_b2 = 0.2221 #[m]
        
        #Eccentric force distance
        l_ec = 0.399 #[m]
        
        #Mass crank
        rho_crank = 8000 #kg/m3
        l_c = 0.200 #[m]
        m_ec = rho_crank*(pi/2*(Re-Rs)**2*l_c)
        
        F_ec = (m_r+m_ec)*e*omega**2
        F_b2 = F_ec*(l_ec-l_b1)/(l_b1-l_b2)
        F_b1 = F_ec + F_b2
        
        #Estimation of balancers
        e_b1 = 0.03 #[m]
        m_b1 = e/e_b1*(m_r + m_ec)*(1+ (l_b1-l_ec)/(l_b2-l_b1))
        e_b2 = 0.03 #[m]
        m_b2 = e/e_b2*(m_r + m_ec)*(l_b1-l_ec)/(l_b2-l_b1)
        
    def friction_losses(self,shell_pressure = 'mean:shell'):

        """
        Friction losses are caused by rubbing parts, bearings etc
        
        Ooi K.T. 'Design optimization of a rolling piston compressor 
        for refrigeration', Applied Therm. Eng. 25(2005) 813-829

        Parameters
        ----------
            shell_pressure : string, 'low', 'low:shell', 'mid', or 'high'
        
            low uses the upstream pressure of the machine,
            
            low:shell uses the pressure after the inlet tube
            
            mid uses the average of upstream and downstream pressures
            
            high uses the pressure downstream of the machine 

        """
        #inlet pressure [Pa]
        inlet_pressure = self.Tubes.Nodes[self.key_inlet].p*1000
        outlet_pressure = self.Tubes.Nodes[self.key_outlet].p*1000
        
        # Find the tube with the inlet node
        Tube = self.Tubes[self.key_inlet]
        # Get the state that is not the inlet state
        if Tube.key1 == 'self.key_inlet':
            shell_pressure_val = Tube.State2.p
        else:
            shell_pressure_val = Tube.State1.p
        
        # Get the shell pressure based on either the inlet or outlet pressure
        # based on whether it is a low-pressure or high-pressure shell
        if shell_pressure == 'low':
            back_pressure = min((inlet_pressure, outlet_pressure))
        elif shell_pressure == 'low:shell':
            back_pressure = min((shell_pressure_val, outlet_pressure))
        elif shell_pressure == 'high':
            back_pressure = max((inlet_pressure, outlet_pressure))
        elif shell_pressure == 'mid':
            back_pressure = (inlet_pressure + outlet_pressure)/2
        elif shell_pressure == 'mean:shell':         
        #Average Pressure shell top of vane
            back_pressure = np.mean(self.CVs['shell'].State.p*1000)
        else:
            raise KeyError("keyword argument shell_pressure must be one of 'low', 'low:shell', 'mid' or 'high'; received '"+str(shell_pressure)+"'")

        self.calculate_force_terms(p_shell = back_pressure)
        self.Wdot_friction = self.Ler + self.Lrc + self.Lec + self.Lv + self.Ls + self.Lbs + self.Lbl
        
        Nstep = self.Itheta+1
        _slice = range(self.Itheta+1)
        
        self.Ler_mean = (np.trapz(self.Ler, self.t[_slice])/(2*pi))
        self.Lrc_mean = (np.trapz(self.Lrc, self.t[_slice])/(2*pi))
        self.Lec_mean =(np.trapz(self.Lec, self.t[_slice])/(2*pi))
        self.Lv_mean = (np.trapz(self.Lv, self.t[_slice])/(2*pi))
        self.Ls_mean = (np.trapz(self.Ls, self.t[_slice])/(2*pi))
        self.Lbs_mean = (np.trapz(self.Lbs, self.t[_slice])/(2*pi))
        self.Lbl_mean = (np.trapz(self.Lbl, self.t[_slice])/(2*pi))
        
        self.Wdot_friction_mean = (np.trapz(self.Wdot_friction, self.t[_slice])/(2*pi))/1000
        print('average mechanical and friction losses [kW]:',self.Wdot_friction_mean)
        
        return self.Wdot_friction_mean

    def mechanical_losses(self):

        if hasattr(self.mech,'detailed_analysis') and self.mech.detailed_analysis == True:            
            #Friction losses
            Wdot_friction = self.friction_losses() #[kW]
            return Wdot_friction
        else:
            print('Simplified mechanical losses')
            return self.mech.Wdot_parasitic #[kW] 
    
    def ambient_heat_transfer(self):
        """
        The ambient heat transfer for the compressor in kW
        
        Returns a positive value if heat is added to the compressor from the 
        ambient
        
        The heat transfer is obtained for natural convection from a vertical cylinder
        """            
        ref = 'Air'
        T_w = self.Tlumps[0] # [K]
        T = self.Tamb # [K]
        P_film = 101.325 #[kPa]
        
        # Vertical side shell
        h_c_vertical = HTC('vertical_plate',T_w,T,P_film,ref,self.Hshell)
        A_shell_vertical = 2*pi*self.Rshell*self.Hshell # m2
        Q_vertical = h_c_vertical*A_shell_vertical*(T-T_w) 

        # Top Surface
        h_c_topshell = HTC('horizontal_plate',T_w,T,P_film,ref,self.Rshell*2,PlateNum='upper_surface')
        A_shell_top = pi*self.Rshell**2 # m2
        Q_topshell = h_c_topshell*A_shell_top*(T-T_w)                     

        # Bottom Surface
        h_c_bottomshell = HTC('horizontal_plate',T_w,T,P_film,ref,self.Rshell*2,PlateNum='lower_surface')
        A_shell_bottom = pi*self.Rshell**2 # m2
        Q_bottomshell = h_c_bottomshell*A_shell_bottom*(T-T_w)  
        
        # Total heat transfer
        Q = Q_vertical + Q_topshell + Q_bottomshell # [kW]
        
        return Q
        
    def lump_energy_balance_callback(self):
        """
        A callback used in PDSimCore.solve to do the energy balance on the lump
        
        Note: we neglect heat transfer to the gas in the working chamber
        """
        #Mechanical and friction losses are added to the lump
        self.Wdot_losses = self.mechanical_losses() #[kW]
        
        #Heat transfer tubes with shell
        self.Qtubes = -sum([Tube.Q for Tube in self.Tubes])
        
        #Heat transfer between the shell and the ambient
        self.Qamb = self.ambient_heat_transfer() #[kW]
        
        #Mechanical power
        self.Wdot_mechanical = self.Wdot_pv + self.Wdot_losses
        
        #The actual torque required to do the compression [N-m]
        self.tau_mechanical = self.Wdot_mechanical / self.omega * 1000
        
        #Mechanical efficiency
        self.eta_mechanical = self.Wdot_pv/self.Wdot_mechanical
        
        if self.motor.type == 'const_eta_motor':
            self.eta_motor = self.motor.eta_motor
        elif self.motor.type == 'motor_map':
            # Use the motor map to calculate the slip rotational speed [rad/s]
            # and the motor efficiency as a function of the torque [N-m]
            eta, omega = self.motor.apply_map(self.tau_mechanical)
            self.eta_motor = eta
            self.omega = omega
        else:
            raise AttributeError('motor.type must be one of "const_eta_motor" or "motor_map"')
        
        #Motor losses [kW]
        self.motor.losses = self.Wdot_mechanical*(1/self.eta_motor-1)
        
        #Electrical Power
        self.Wdot_electrical = self.Wdot_mechanical + self.motor.losses
        
        if hasattr(self,'Wdot_i'):
            #Overall isentropic efficiency
            self.eta_oi = self.Wdot_i/self.Wdot_electrical        

        if self.verbosity > 0:
            print('At this iteration')
            print('    Indicated power:', self.Wdot_pv,'kW')
            print('    Mechanical torque:',self.tau_mechanical,'Nm')
            print('    Electrical power:', self.Wdot_electrical,'kW')
            print('    Mechanical losses:', self.Wdot_losses,'kW')
            print('    Motor losses:', self.motor.losses,'kW')
            print('    Mass flow rate:', self.mdot,'kg/s', self.mdot*3600,'kg/h')
            if hasattr(self,'Wdot_i'):
                print('    Over. isentropic:', self.eta_oi,'-')
            if hasattr(self,'eta_v'):
                print('    Volumetric:', self.eta_v,'-')
            if self.mech.detailed_analysis:
                print('Mechanical efficiency:',self.eta_mechanical,'-')
                print('Ecc. & inner roller :', self.Ler_mean,'W')
                print('Roller face & cylinder head :', self.Lrc_mean,'W')
                print('Ecc. face & cylinder head :', self.Lec_mean,'W')
                print('Vane tip & roller :', self.Lv_mean,'W')
                print('Vane sides & vane slot :', self.Ls_mean,'W')
                print('Shaft & upper bearing :', self.Lbs_mean,'W')
                print('Shaft & lower bearing :', self.Lbl_mean,'W')
            
        return  self.Qamb + self.Qtubes + 0.1*self.motor.losses

    def step_callback(self,t,h,Itheta):
        """
        Here we test whether the control volumes need to be
        a) Merged
        b) Adjusted because you are at the discharge angle
        
        """ 
        if hasattr(self.mech,'cylinder_config') and self.mech.cylinder_config == 'dual-cylinder':
            if self.solver == 'RK45':
            
                # This gets called at every step, or partial step
                self.theta = t
                
                def angle_difference(angle1,angle2):
                    # Due to the periodicity of angles, you need to handle the case where the
                    # angles wrap around - suppose theta_d is 6.28 and you are at an angles of 0.1 rad
                    #, the difference should be around 0.1, not -6.27
                    # 
                    # This brilliant method is from http://blog.lexique-du-net.com/index.php?post/Calculate-the-real-
                    #difference-between-two-angles-keeping-the-sign
                    # and the comment of user tk
                    return (angle1-angle2+pi)%(2*pi)-pi
                
                disable=False
                
                if t<pi<t+h: #and t<pi - 2e-15:
                    # Take a step almost up to the discharge angle
                    # print('almost at discontinuity')
                    disable = True
                    h = pi-t-1e-15                    
                elif abs(angle_difference(pi,t)) < 1e-14:                    
                    for key,angle in self.discontinuity_angle.items():
                        #oldCVs
                        oldCV_C=self.CVs['C']
                        oldCV_D=self.CVs['D']
                        # 
                        if key == 'C':
                            #Set the state of the "new" chamber to be the old chamber
                            self.CVs[key].State.update({'T':oldCV_D.State.T,'D':oldCV_D.State.rho})
    
                        elif key == 'D':
                            #Set the state of the "new" chamber to be the old chamber
                            self.CVs[key].State.update({'T':self.Tin,'D':self.rhoin})
    
                    #Re-calculate the CV volumes
                    V,dV = self.CVs.volumes(pi+1e-10)
        
                    #Update the matrices using the new CV definitions
                    self.T[self.CVs.exists_indices,Itheta] = self.CVs.T
                    self.p[self.CVs.exists_indices,Itheta] = self.CVs.p
                    self.m[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)*V
                    self.rho[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)
                    
                    #print V
                    # This has to be quite large - why?
                    h = 2e-10
                    disable='no_integrate'
                elif t > pi:                    
                    #print 'here h:',h,t
                    disable = False                
        #         if t > 3.5:
        #             raise ValueError                
                return disable,h    
            else:
                self.theta = t
                return False, h
        else:
            self.theta = t
            return False, h


    

