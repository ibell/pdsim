## Python imports
from math import pi
import numpy as np
from scipy.optimize import newton
from matplotlib import pyplot as plt
import matplotlib.pyplot as plt


## PDSim imports
from PDSim.core.core import PDSimCore
from PDSim.core.containers import Tube,ControlVolume
from PDSim.core.motor import Motor

import PDSim.scroll.core as core
from PDSim.scroll import scroll_geo
from PDSim.scroll.plots import plotScrollSet

from PDSim.flow import flow_models
from PDSim.flow.flow import FlowPath
from PDSim.flow.flow_models import IsentropicNozzleWrapper

from PDSim.misc.datatypes import arraym

from PDSim.plot.plots import debug_plots

from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import State as CPState
from PDSim.core.state_flooded import StateFlooded


class ScrollExpander(core.Scroll):
    
    #The expander is just like the compressor, except that the volumes go backwards. 
    #Theta starts a revolution at 2*pi and then decreases over one rotation to 0. 
    #Beta is equal to 2*pi-theta whc.
    
    
    
    def __init__(self, *args, **kwargs):
        core.Scroll.__init__(self, *args, **kwargs)
        self.__before_discharge1__ =True 
        self.__before_discharge2__ =True
        self.__hasLiquid__ = False
        
    def V_ss(self, beta, full_output=False):
        VdV = core.Scroll.V_dd(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
        
    def V_s1(self, beta, full_output=False):
        VdV = core.Scroll.V_d1(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
    
    def V_s2(self, beta, full_output=False):
        VdV = core.Scroll.V_d2(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
    
    def V_sss(self, beta, full_output=False):
        if self.__before_discharge1__ == True and beta < 2*pi - self.theta_d:
            VdV = core.Scroll.V_dd(self, 2*pi-beta, full_output=full_output)
        else:
            VdV = core.Scroll.V_ddd(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
        
    def V_e1(self, beta, alpha=1, full_output=False):
        VdV = core.Scroll.V_c1(self, 2*pi-beta, alpha=alpha, full_output=full_output)
        return VdV[0],-VdV[1]
    
    def V_e2(self, beta, alpha=1, full_output=False):
        VdV = core.Scroll.V_c2(self, 2*pi-beta, alpha=alpha, full_output=full_output)
        return VdV[0],-VdV[1]
    
    def V_d1(self, beta):
        VdV = core.Scroll.V_s1(self, 2*pi-beta)
        return VdV[0],-VdV[1]
    
    def V_d2(self, beta):
        VdV = core.Scroll.V_s2(self, 2*pi-beta)
        return VdV[0],-VdV[1]
    
    def V_da(self, beta, full_output=False):
        VdV = core.Scroll.V_sa(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
    
## Adds all the control volumes for the scroll expander 
    def auto_add_CVs(self,inletState,outletState):  
    
     
        #``s1``        Suction chamber on side 1
        #``s2``        Suction chamber on side 2
        #``ss``        Central suction chamber
        #``sss``       Merged suction chamber
        #``e1.i``      The i-th expansion chamber on side 1 (i=1 for outermost chamber)
        #``e2.i``      The i-th expansion chamber on side 2 (i=1 for outermost chamber)
        #``da``        Discharge Area
        #``d1``        Discharge chamber on side 1
        #``d2``        Discharge chamber on side 2
    
        #Suction & Expansion chambers
        self.add_CV(ControlVolume(key = 'sss', initialState = inletState.copy(),VdVFcn = self.V_sss))
                
        self.add_CV(ControlVolume(key = 'ss', initialState = inletState.copy(),VdVFcn = self.V_ss, exists = False, discharge_becomes = 'sss'))
        
        self.add_CV(ControlVolume(key = 's1', initialState = inletState.copy(),VdVFcn = self.V_s1, exists = False, discharge_becomes = 'e1.'+str(scroll_geo.nC_Max(self.geo))))
        
        self.add_CV(ControlVolume(key = 's2', initialState = inletState.copy(),VdVFcn = self.V_s2, exists = False, discharge_becomes = 'e2.'+str(scroll_geo.nC_Max(self.geo))))
        
        #Discharge chambers
        self.add_CV(ControlVolume(key='da', initialState = outletState.copy(),VdVFcn=self.V_da))
        
        self.add_CV(ControlVolume(key='d1', initialState = outletState.copy(),VdVFcn=self.V_d1,becomes='none'))
        
        self.add_CV(ControlVolume(key='d2', initialState = outletState.copy(),VdVFcn=self.V_d2,becomes='none'))
        
        #Add each pair of expansion chambers
        nCmax = scroll_geo.nC_Max(self.geo)
        
        #Must have at least one pair
        assert (nCmax>=1)
        
        #Innermost pair at nCmax - 1 starts at inlet state
        for alpha in range(nCmax,0,-1): 
            keye1 = 'e1.'+str(alpha)
            keye2 = 'e2.'+str(alpha)
            if alpha == nCmax:
                #It is the outermost pair of expansion chambers
                if self.__hasLiquid__ == False:
                    initState = CPState(inletState.Fluid,dict(T=inletState.T,D=inletState.rho))
                elif self.__hasLiquid__ == True: 
                    initState = StateFlooded(inletState.Fluid,inletState.Liq,inletState.p,inletState.T,self.xL,'HEM')
                    initState.update(dict(T=inletState.T, P = inletState.p, xL=self.xL))
                else:
                    raise Exception('Not implemented')
                    
            else:
                #It is not the first CV, more involved analysis
                #Assume isentropic expansion from the inlet state at the end of the suction process
                T1 = inletState.T
                p1 = inletState.p
                s1 = inletState.s
                rho1 = inletState.rho
                k = inletState.cp/inletState.cv#Specific heat ratio
                V1 = self.V_sss(2*pi-self.theta_d-1e-14)[0]
                V2 = self.V_e1(0, alpha)[0]*2
                # Mass is constant, so rho1*V1 = rho2*V2
                rho2 = rho1 * V1 / V2
                # Now don't know temperature or pressure, but you can assume it is isentropic to find the temperature
                if self.__hasLiquid__ == False:
                    T2 = newton(lambda T: PropsSI('S','T',T,'D',rho2,inletState.Fluid)/1000-s1, T1)
                    initState = CPState(inletState.Fluid,dict(T=T2,D=rho2)).copy()
                elif self.__hasLiquid__ == True: 
                    #TODO: add a better guess
                    kstar_mix = inletState.get_kstar()
                    T2s = T1*(outletState.p/p1)**((kstar_mix - 1)/kstar_mix)
                    T2 = T1 - (T1-T2s)*0.7
                    initState = StateFlooded(inletState.Fluid,inletState.Liq,outletState.p,T2,inletState.xL,'HEM').copy()
                    #(initState.update(dict(T=T2, P = outletState.p, xL=self.xL))).copy()

                else:
                    raise Exception('Not implemented')                
            
            # Expansion chambers do not change definition at discharge angle
            disc_becomes_e1 = keye1
            disc_becomes_e2 = keye2
            
            if alpha > 1:
                exists = False # Starts out not being in existence, comes into existence at the discharge angle
            else:
                exists = True 
                               
            if alpha == nCmax:
                # It is the innermost pair of chambers, becomes a discharge chamber at the end of the rotation
                becomes_e1 = 'e1.'+str(alpha-1)
                becomes_e2 = 'e2.'+str(alpha-1)
            else:
                # It is not the innermost pair of chambers, becomes another set of expansion chambers at the end of the rotation
                becomes_e1 = 'd1'
                becomes_e2 = 'd2'
                               
            self.add_CV(ControlVolume(key = keye1,initialState = initState.copy(),VdVFcn = self.V_e1,VdVFcn_kwargs = {'alpha':alpha},discharge_becomes = disc_becomes_e1,becomes = becomes_e1,exists = exists))
            
            self.add_CV(ControlVolume(key = keye2,initialState = initState.copy(),VdVFcn = self.V_e2,VdVFcn_kwargs = {'alpha':alpha},discharge_becomes = disc_becomes_e2,becomes = becomes_e2,exists = exists))
     
## Callback - Energy balance on the lumps       
    def lump_energy_balance_callback(self):
                
        self.Qamb = self.ambient_heat_transfer(self.Tlumps[0])
        
        self.Qsuc = sum([Tube.Q for Tube in self.Tubes])
        
        self.Qscroll = self.HTProcessed.mean_Q
        
        #self.mech.Wdot_losses = self.mechanical_losses('low:shell') 
        #self.mech.Wdot_losses = (1*2*3.14*38.25)/1000.0
        self.mech.Wdot_losses = 0.5*self.omega/1000.0
        
        #Shaft power from forces on the orbiting scroll from the gas in the pockets [kW]
        #self.Wdot_forces = self.omega*self.forces.mean_tau

        #HT terms are positive if heat transfer is to the lump
        Qnet = 0.0   
        
        #Suction heating
        Qnet -= self.Qsuc

        #Heat transfer with the ambient
        #Qnet -= self.Qamb   #(Qamb is positive if heat is being removed, so flip the)
        Qnet += self.Qamb
        
        #Mechanical losses
        Qnet -= self.mech.Wdot_losses
        
        # Heat transfer with the gas in the working chambers 
        #Qnet -= self.Qscroll #(mean_Q is positive if heat is transfered to the gas in the working chamber, so flip the sign for the lump)
        Qnet += self.Qscroll
        
        #self.Res = self.Qamb + self.Qscroll - self.mech.Wdot_losses - self.Qsuc
        
        #Total shaft work
        self.h1 = inletState.h
        self.h2 = outletState.h
        self.Wdot = self.mdot*(self.h1-self.h2) - self.Qamb
        
        self.Wdot_shaft = self.Wdot - self.mech.Wdot_losses
            
        #Total mechanical work
        self.Wdot_mechanical = self.Wdot_pv + self.mech.Wdot_losses
        
        #The actual torque generated by expansion [N-m]
        self.tau_mechanical = self.Wdot_mechanical / self.omega * 1000  
        
        # Motor losses - Constant efficiency or Motor map
        if self.motor.type == 'const_eta_motor':
            self.eta_motor = self.motor.eta_motor
        elif self.motor.type == 'motor_map':    
            eta, omega = self.motor.apply_map(self.tau_mechanical)
            self.eta_motor = eta
            self.omega = omega
        else:
            raise AttributeError('motor.type must be one of "const_eta_motor" or "motor_map"')
        
        #Motor losses [kW]
        #self.motor.losses = self.Wdot_mechanical*(1/self.eta_motor-1) 
        self.motor.losses = self.Wdot_shaft*(1 - self.eta_motor)
        
        #Electrical Power
        #self.Wdot_electrical = self.Wdot_mechanical + self.motor.losses
        self.Wdot_electrical = self.Wdot_shaft - self.motor.losses 
        
        if hasattr(self,'Wdot_i'):
            self.eta_oi = self.Wdot_electrical/self.Wdot_i  #Overall isentropic efficiency
        
        print 'At this iteration'
        print 'Electrical power:', self.Wdot_electrical,'kW'
        print 'Mass flow rate:', self.mdot,'kg/s'
        
        if hasattr(self,'Wdot_i'):
            print 'Over. isentropic:', self.eta_oi,'-'
            
        if hasattr(self,'eta_v'):
            print 'Volumetric:', self.eta_v,'-'
            
        return [Qnet]
        
 ## Flow paths   
    def S_to_SS(self, FlowPath, **kwargs):
        if self.__before_discharge1__:
            FlowPath.A = 0.0
        else:
            FlowPath.A = scroll_geo.Area_d_dd(2*pi-self.beta, self.geo)
        try:
            mdot = flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
        
    def DA_D(self, FlowPath, X_d=1.0,**kwargs):
        FlowPath.A=X_d*scroll_geo.Area_s_sa(2*pi-self.beta, self.geo)
        try:
            return flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
        except ZeroDivisionError:
            return 0.0
    
## Test - whether control volumes needs to be a) Split b) Adjusted at the discharge angle 
       
    def step_callback(self,t,h,Itheta):
        #This gets called at every step, or partial step
        self.beta=t
        
        def angle_difference(angle1,angle2):
        #Due to the periodicity of angles, you need to handle the case where the angles wrap around
        #Suppose theta_d is 6.28 and you are at an angles of 0.1 rad, the difference should be around 0.1, not -6.27. 
            
            return (angle1-angle2+pi)%(2*pi)-pi
        
        def IsAtSplit():
            if t > 3.5 and t < 2*pi-self.theta_d:
                return True
            else:
                return False
            
        disable=False
        
        if t<2*pi-self.theta_d<t+h and self.__before_discharge2__==False:
            #Take a step almost up to the discharge angle using the old definitions
            disable=True
            h=2*pi-self.theta_d-t-1e-10
            self.__before_discharge2__=True
        elif self.__before_discharge2__==True:
            
            #At the discharge angle
            print 'At the discharge angle'
            
            #Reassign chambers
            #Find chambers with a discharge_becomes flag
            for key in self.CVs.exists_keys:
                if self.CVs[key].discharge_becomes in self.CVs.keys:
                    #Set the state of the "new" chamber to be the old chamber
                    oldCV=self.CVs[key]
                    #print 'setting',key,'to',oldCV.discharge_becomes
                    if oldCV.exists==True:
                        newCV=self.CVs[oldCV.discharge_becomes]
                        newCV.State.update({'T':oldCV.State.T,'D':oldCV.State.rho})
                        oldCV.exists = False
                        newCV.exists = True
                    else:
                        raise AttributeError("old CV doesn't exist")
            
            self.__before_discharge2__ = False
            self.__before_discharge1__ = True
            
            self.update_existence()
            
            #Re-calculate the CV volumes
            V,dV = self.CVs.volumes(t)
            #Update the matrices using the new CV definitions
            self.T[self.CVs.exists_indices,Itheta]=self.CVs.T
            self.p[self.CVs.exists_indices,Itheta]=self.CVs.p
            self.m[self.CVs.exists_indices,Itheta]=arraym(self.CVs.rho)*V
            self.rho[self.CVs.exists_indices,Itheta]=arraym(self.CVs.rho)
            
            # Adaptive makes steps of h/4 3h/8 12h/13 and h/2 and h
            # Make sure step does not hit any *right* at theta_d
            # That is why it is 2.2e-10 rather than 2.0e-10
            h=2.2e-10
            disable=True
       
        elif self.CVs['sss'].exists and IsAtSplit():
            
            #Build the volume vector using the old set of control volumes (pre-merge)
            V,dV=self.CVs.volumes(t)
            
            print 'splitting'
            
            if self.__hasLiquid__== False:

                T = self.CVs['sss'].State.T
                rho = self.CVs['sss'].State.rho
                self.CVs['s1'].State.update({'T':T,'D':rho})
                self.CVs['s2'].State.update({'T':T,'D':rho})
                self.CVs['ss'].State.update({'T':T,'D':rho})
                
                self.CVs['sss'].exists = False
                self.CVs['s1'].exists = True
                self.CVs['s2'].exists = True
                self.CVs['ss'].exists = True
                
                self.update_existence()
                
                #Re-calculate the CV
                V,dV=self.CVs.volumes(t)
                self.T[self.CVs.exists_indices,Itheta] = self.CVs.T
                self.p[self.CVs.exists_indices,Itheta] = self.CVs.p
                self.m[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)*V
                self.rho[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)
                
            else:
                raise NotImplementedError('no flooding yet')
            
            h = 1e-8
            disable=True 
              
        elif t > 2*pi - self.theta_d:
            self.__before_discharge1__ = False
            disable = False
            
        return disable, h
    
    def pre_run(self):
        """
        Intercepts the call to pre_run and does some scroll processing, then 
        calls the base class function
        """
        
        #Call the base class function        
        PDSimCore.pre_run(self)
                
#def Expander(Tinlet,pinlet,poutlet,N,Tamb):
def Expander():    
    SE = ScrollExpander()
    
## Model Inputs
    
    Ref = 'R123'    #Refrigerant
    Tinlet = 394    #Inlet temperature
    pinlet = 545    #Inlet pressure
    poutlet = 200   #Outlet pressure
    N = 2295
    Tamb = 304
    Liq = 'Duratherm_LT'
    SE.xL = 0.5   #Liquid fraction
    
    #Flood_Prop = State_Flooded(Ref,Liq,Tinlet,pinlet,xL)
    #SE.Entropy = Flood_Prop.s_mix(Ref,Liq,Tinlet,pinlet,xL)
    
    #SE.h_mix_in = h_m(Ref,"POE",Tinlet,pinlet,xL)
    #SE.rho_mix_in = rho_m(Ref,"POE",Tinlet,pinlet,xL)
    #SE.v_mix_in = 1/SE.rho_mix_in
    #Guess outlet temperature
    #SE.Toutlet_s = newton(lambda T:s_m(Ref,"POE",T,poutlet,xL)-SE.s_mix_in,Tinlet)
    #Mixture properties at inlet
    #SE.s_mix_out = s_m(Ref,"POE",SE.Toutlet_s,poutlet,xL)
    #SE.h_mix_out = h_m(Ref,"POE",SE.Toutlet_s,poutlet,xL)
    #SE.rho_mix_out = rho_m(Ref,"POE",SE.Toutlet_s,poutlet,xL)
    #SE.v_mix_out = 1/SE.rho_mix_out
    #SE.mdot_guess_mix = SE.rho_mix_in*SE.Vdisp/SE.Vratio*SE.omega/(2*pi)
    
    SE.Tamb = Tamb #Ambient temperature
    SE.omega = (N*2*pi)/60 #Rotational speed

    SE.set_scroll_geo(110e-6,2.24,0.0046,0.0064)   #scroll geometry
    SE.set_disc_geo('2Arc', r2 = 0) #Discharge geometry
    SE.geo.delta_radial = 0
    SE.geo.delta_flank = 70e-6
    
    SE.geo.delta_suction_offset = 0.0e-3
    SE.geo.phi_ie_offset = 0.0

    SE.mech = core.struct() #Set the bearing dimensions
    SE.mech.D_upper_bearing = 0.02540
    SE.mech.L_upper_bearing = 0.031
    SE.mech.c_upper_bearing = 10e-6
    SE.mech.D_crank_bearing = 0.0254
    SE.mech.L_crank_bearing = 0.023
    SE.mech.c_crank_bearing = 10e-6
    SE.mech.D_lower_bearing = 0.01895
    SE.mech.L_lower_bearing = 0.023
    SE.mech.c_lower_bearing = 10e-6 
    
    SE.mech.thrust_ID = 0.05    #Thrust bearing inner diameter
    SE.mech.thrust_friction_coefficient = 0.02 #Tuned for point A1
    SE.mech.L_ratio_bearings = 5    #Ratio of the distances from the upper bearing to the crank bearing
    SE.mech.mu_oil = 0.008  #Viscosity of the oil
    
    SE.mech.scroll_plate_thickness = 0.007 #Scroll plate thickness
    SE.mech.scroll_density = 8100   #Scroll density
    SE.mech.scroll_added_mass = 0   #Added mass
    SE.mech.scroll_plate_diameter = 0.092 #Scroll plate diameter
    m, zcm = SE.calculate_scroll_mass()
    SE.mech.orbiting_scroll_mass = m    #Scroll mass
    SE.mech.scroll_zcm__thrust_surface = zcm    #Centroid location
    
    SE.h_shell = 0.02   #Shell heat transfer coefficient
    SE.A_shell = 0.05   #Shell area for ambient heat transfer
    
    SE.HTC = 1
    
    SE.motor = Motor()
    SE.motor.set_eta(0.9)   #Motor efficiency
    SE.motor.suction_fraction = 0.0 # (If fraction of heat from motor losses get added to suction flow, by adding Tube.Q_add to Tube.Q)
    
    SE.mech.detailed_analysis = True
    
## Inlet/Outlet State
    global inletState
    global outletState


    if SE.__hasLiquid__ == False:
    
        inletState = CPState(Ref,{'T':Tinlet,'P':pinlet})
        T2s = SE.guess_outlet_temp(inletState,poutlet)  #Guess outlet
        outletState = CPState(Ref,{'T':T2s,'P':poutlet})
    elif SE.__hasLiquid__ == True:
        inletState = StateFlooded(Ref,Liq,pinlet,Tinlet,SE.xL,'HEM')
        inletState.update(dict(T=Tinlet, P=pinlet, xL=SE.xL))
        T2s = SE.guess_outlet_temp(inletState,poutlet)  #Guess outlet
        outletState = StateFlooded(Ref,Liq,poutlet,T2s,SE.xL,'HEM')        
        outletState.update(dict(T=T2s, P=poutlet, xL=SE.xL))
        
## Define CVc and tubes - Suction, Expansion and Exhaust
    #Guess mass flow rate
    SE.mdot_guess = inletState.rho*SE.Vdisp/SE.Vratio*SE.omega/(2*pi)
    
    
    if SE.__hasLiquid__ == False:
        #Suction tube
        SE.add_tube(Tube(key1 = 'inlet.1',key2 = 'inlet.2',L = 0.086,ID = 0.018,mdot = SE.mdot_guess, State1 = inletState.copy(),fixed = 1,TubeFcn = SE.TubeCode))
        #Exhaust tube
        SE.add_tube(Tube(key1 = 'outlet.1',key2 = 'outlet.2',L = 0.09,ID = 0.015,mdot = SE.mdot_guess, State2 = outletState.copy(),fixed = 2,TubeFcn = SE.TubeCode))
    
    elif SE.__hasLiquid__ == True:
        #Suction tube
        SE.add_tube(Tube(key1 = 'inlet.1',key2 = 'inlet.2',L = 0.086,ID = 0.018,mdot = SE.mdot_guess, StateFlood1 = inletState.copy(),fixed = 1,TubeFcn = SE.TubeCode,__hasLiquid__ = SE.__hasLiquid__))
        #Exhaust tube
        SE.add_tube(Tube(key1 = 'outlet.1',key2 = 'outlet.2',L = 0.09,ID = 0.015,mdot = SE.mdot_guess, StateFlood2 = outletState.copy(),fixed = 2,TubeFcn = SE.TubeCode,__hasLiquid__ = SE.__hasLiquid__))        
    else:
        print 'Not implemented'
    
    #Adds all the control volumes for the scroll expander
    SE.auto_add_CVs(inletState, outletState)
    
## Add all the leakage terms for the expander
    SE.auto_add_leakage(flankFunc = SE.FlankLeakage,radialFunc = SE.RadialLeakage)
    
## Add flows between CVs and tubes

    FP = FlowPath(key1='inlet.2',key2='ss',MdotFcn=IsentropicNozzleWrapper(),)
    FP.A = pi*0.007**2/4
    #FP.A = pi*0.018**2/4
    SE.add_flow(FP)
    
    FP = FlowPath(key1='inlet.2',key2='sss',MdotFcn=IsentropicNozzleWrapper(),)
    FP.A = pi*0.007**2/4
    #FP.A = pi*0.018**2/4
    SE.add_flow(FP)
    
    SE.add_flow(FlowPath(key1='s1',key2='ss',MdotFcn=SE.S_to_SS))
    
    SE.add_flow(FlowPath(key1='s2',key2='ss',MdotFcn=SE.S_to_SS))

    FP = FlowPath(key1='outlet.1',key2='da',MdotFcn=IsentropicNozzleWrapper(),)
    FP.A = pi*0.02**2/4
    #FP.A = pi*0.015**2/4
    SE.add_flow(FP)
    
    SE.add_flow(FlowPath(key1 = 'da',key2 = 'd1',MdotFcn = SE.DA_D,MdotFcn_kwargs = dict(X_d = 0.5)))
    
    SE.add_flow(FlowPath(key1 = 'da',key2 = 'd2',MdotFcn = SE.DA_D,MdotFcn_kwargs = dict(X_d = 0.5)))
    
    SE.Qwall = SE.wrap_heat_transfer
    
## Calculations
    SE.connect_callbacks(step_callback = SE.step_callback, heat_transfer_callback = SE.heat_transfer_callback,lumps_energy_balance_callback = SE.lump_energy_balance_callback,endcycle_callback= SE.endcycle_callback)

    SE.RK45_eps = 1e-6
    SE.eps_cycle = 0.003
    
    SE.precond_solve(key_inlet='inlet.1',key_outlet='outlet.2',solver_method='RK45',UseNR = False,OneCycle = False,plot_every_cycle= False,hmin = 1e-8)
    
    #del SE.FlowStorage
    #from PDSim.misc.hdf5 import HDF5Writer
    #h5 = HDF5Writer()
    #h5.write_to_file(SE, 'exp_data.h5')

    mass = SE.mdot
    work = SE.Wdot_shaft
    T_ex = SE.T2
    efficinecy = SE.Wdot_shaft/SE.Wdot_i
    #Test_value = SE.Entropy

    return [mass,work,T_ex,efficinecy]
     
if __name__=='__main__':
    Expander()
