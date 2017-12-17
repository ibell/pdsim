from __future__ import print_function

from math import pi
import numpy as np
from PDSim.scroll import scroll_geo
import PDSim.scroll.core as core
from PDSim.core.motor import Motor
from PDSim.flow.flow import FlowPath
from PDSim.flow.flow_models import IsentropicNozzleWrapper
from PDSim.core.containers import Tube,ControlVolume
import matplotlib.pyplot as plt
from CoolProp.CoolProp import Props
from CoolProp.State import State as CPState
from PDSim.core.core import PDSimCore
from PDSim.flow import flow_models
from PDSim.misc.datatypes import arraym

# If scipy is available, use its interpolation and optimization functions, otherwise, 
# use our implementation (for packaging purposes mostly)
try:
    import scipy.optimize as optimize
except ImportError:
    import PDSim.misc.solvers as optimize

class ScrollExpander(core.Scroll):
    """
    The expander is just like the compressor, except that the volumes go backwards.
    
    theta starts a revolution at 2*pi, and then decreases over one rotation to 0. Beta is equal to 2*pi-theta whc
    In this way, all the code for the compressor can be more readily used with the expander
    """
    
    def __init__(self, *args, **kwargs):
        core.Scroll.__init__(self, *args, **kwargs)
        
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
    
    def V_s1(self, beta, full_output=False):
        VdV = core.Scroll.V_d1(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
    
    def V_s2(self, beta, full_output=False):
        VdV = core.Scroll.V_d2(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
    
    def V_ss(self, beta, full_output=False):
        VdV = core.Scroll.V_dd(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
    
    def V_sss(self, beta, full_output=False):
        if self.__before_discharge1__ == True and beta < 2*pi - self.theta_d:
            VdV = core.Scroll.V_dd(self, 2*pi-beta, full_output=full_output)
        else:
            VdV = core.Scroll.V_ddd(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
    
    def V_da(self, beta, full_output=False):
        VdV = core.Scroll.V_sa(self, 2*pi-beta, full_output=full_output)
        return VdV[0],-VdV[1]
    
    def auto_add_CVs(self,inletState,outletState):
        """
        Adds all the control volumes for the scroll expander.
        
        Parameters
        ----------
        inletState
            A :class:`State <CoolProp.State.State>` instance for the inlet to the scroll set.  Can be approximate
        outletState
            A :class:`State <CoolProp.State.State>` instance for the outlet to the scroll set.  Can be approximate
            
        Notes
        -----
        Uses the indices of 
        
        ============= ===================================================================
        CV            Description
        ============= ===================================================================
        ``da``        Discharge Area
        ``d1``        Discharge chamber on side 1
        ``d2``        Discharge chamber on side 2
        ``s1``        Suction chamber on side 1
        ``s2``        Suction chamber on side 2
        ``ss``        Central suction chamber
        ``sss``       Merged suction chamber
        ``e1.i``      The i-th expansion chamber on side 1 (i=1 for outermost chamber)
        ``e2.i``      The i-th expansion chamber on side 2 (i=1 for outermost chamber)
        ============= ===================================================================
        """
        
        #Add all the control volumes that are easy.  Suction area and suction chambers
        self.add_CV(ControlVolume(key='da', initialState = outletState.copy(),
                                VdVFcn=self.V_da))
        self.add_CV(ControlVolume(key='d1', initialState = outletState.copy(),
                                VdVFcn=self.V_d1,becomes='none'))
        self.add_CV(ControlVolume(key='d2', initialState = outletState.copy(),
                                VdVFcn=self.V_d2,becomes='none'))
        
        #Discharge chambers are also easy.  Assume that you start with 'ddd' chamber merged.
        # No problem if this isn't true.
        self.add_CV(ControlVolume(key = 's1', initialState = inletState.copy(),
                                  VdVFcn = self.V_s1, exists = False, discharge_becomes = 'e1.'+str(scroll_geo.nC_Max(self.geo))))
        self.add_CV(ControlVolume(key = 's2', initialState = inletState.copy(),
                                  VdVFcn = self.V_s2, exists = False, discharge_becomes = 'e2.'+str(scroll_geo.nC_Max(self.geo))))
        self.add_CV(ControlVolume(key = 'ss', initialState = inletState.copy(),
                                  VdVFcn = self.V_ss, exists = False, discharge_becomes = 'sss'))
        self.add_CV(ControlVolume(key = 'sss', initialState = inletState.copy(),
                                  VdVFcn = self.V_sss))

        #Add each pair of expansion chambers
        nCmax = scroll_geo.nC_Max(self.geo)
        # Must have at least one pair
        assert (nCmax>=1)
        
        # Innermost pair at nCmax - 1 starts at inlet state
        
        for alpha in range(nCmax,0,-1): 
            keye1 = 'e1.'+str(alpha)
            keye2 = 'e2.'+str(alpha)
            if alpha == nCmax:
                #It is the outermost pair of compression chambers
                initState = CPState(inletState.Fluid,
                                        dict(T=inletState.T,
                                             D=inletState.rho)
                                        )
            else:
                #It is not the first CV, more involved analysis
                #Assume isentropic compression from the inlet state at the end of the suction process
                T1 = inletState.T
                s1 = inletState.s
                rho1 = inletState.rho
                k = inletState.cp/inletState.cv
                V1 = self.V_sss(2*pi-self.theta_d-1e-14)[0]
                V2 = self.V_e1(0, alpha)[0]*2
                # Mass is constant, so rho1*V1 = rho2*V2
                rho2 = rho1 * V1 / V2
                # Now don't know temperature or pressure, but you can assume
                # it is isentropic to find the temperature
                T2 = optimize.newton(lambda T: PropsSI('S','T',T,'D',rho2,inletState.Fluid)-s1*1000, T1)
                initState = CPState(inletState.Fluid,dict(T=T2,D=rho2)).copy()
            
            # Expansion chambers do not change definition at discharge angle
            disc_becomes_e1 = keye1
            disc_becomes_e2 = keye2
            
            if alpha > 1:
                exists = False # Starts out not being in existence, comes into 
                               # existence at the discharge angle
            else:
                exists = True 
                               
            if alpha == nCmax:
                # It is the innermost pair of chambers, becomes a discharge
                # chamber at the end of the rotation
                becomes_e1 = 'e1.'+str(alpha-1)
                becomes_e2 = 'e2.'+str(alpha-1)
            else:
                # It is not the innermost pair of chambers, becomes another 
                # set of compression chambers at the end of the rotation
                becomes_e1 = 'd1'
                becomes_e2 = 'd2'
                               
            self.add_CV(ControlVolume(key = keye1,
                                      initialState = initState.copy(),
                                      VdVFcn = self.V_e1,
                                      VdVFcn_kwargs = {'alpha':alpha},
                                      discharge_becomes = disc_becomes_e1,
                                      becomes = becomes_e1,
                                      exists = exists))
            
            self.add_CV(ControlVolume(key = keye2,
                                      initialState = initState.copy(),
                                      VdVFcn = self.V_e2,
                                      VdVFcn_kwargs = {'alpha':alpha},
                                      discharge_becomes = disc_becomes_e2,
                                      becomes = becomes_e2,
                                      exists = exists))
            
    def lump_energy_balance_callback(self):
        """
        The callback where the energy balance is carried out on the lumps
        
        Notes
        -----
        Derivation for electrical power of motor:
        
        .. math ::
            
            \\eta _{motor} = \\frac{\\dot W_{shaft}}{\\dot W_{shaft} + \\dot W_{motor}}
            
        .. math ::
            
            {\\eta _{motor}}\\left( \\dot W_{shaft} + \\dot W_{motor} \\right) = \\dot W_{shaft}
            
        .. math::
        
            \\dot W_{motor} = \\frac{\\dot W_{shaft}}{\\eta _{motor}} - \\dot W_{shaft}
        """
        
        #For the single lump
        # HT terms are positive if heat transfer is TO the lump
        Qnet = 0.0
        Qnet -= sum([Tube.Q for Tube in self.Tubes])
        
        self.Qamb = self.ambient_heat_transfer(self.Tlumps[0])
        self.mech.Wdot_losses = 0.9*self.omega/1000.0
        
        # Heat transfer with the ambient; Qamb is positive if heat is being removed, thus flip the sign
        Qnet -= self.Qamb
        
        Qnet += self.mech.Wdot_losses
        # Heat transfer with the gas in the working chambers.  mean_Q is positive
        # if heat is transfered to the gas in the working chamber, so flip the 
        # sign for the lump
        Qnet -= self.HTProcessed.mean_Q
        
        self.Wdot_mechanical = self.Wdot_pv + self.mech.Wdot_losses
        
        #The actual torque required to do the compression [N-m]
        self.tau_mechanical = self.Wdot_mechanical / self.omega * 1000
        
        # 2 Options for the motor losses:
        # a) Constant efficiency
        # b) Based on the torque-speed-efficiency motor
        
#         if self.motor.type == 'const_eta_motor':
#             self.eta_motor = self.motor.eta_motor
#         elif self.motor.type == 'motor_map':
#             # Use the motor map to calculate the slip rotational speed [rad/s]
#             # and the motor efficiency as a function of the torque [N-m]
#             eta, omega = self.motor.apply_map(self.tau_mechanical)
#             self.eta_motor = eta
#             self.omega = omega
#         else:
#             raise AttributeError('motor.type must be one of "const_eta_motor" or "motor_map"')
        
#         #Motor losses [kW]
#         self.motor.losses = 0#self.Wdot_mechanical*(1/self.eta_motor-1)
        
        #Electrical Power
        self.Wdot_electrical = self.Wdot_mechanical# + self.motor.losses
        
        if hasattr(self,'Wdot_i'):
            #Overall isentropic efficiency
            self.eta_oi = self.Wdot_electrical/self.Wdot_i
        
        print('At this iteration')
        print('    Electrical power:', self.Wdot_electrical,'kW')
        print('    Mass flow rate:', self.mdot,'kg/s')
        if hasattr(self,'Wdot_i'):
            print('    Over. isentropic:', self.eta_oi,'-')
        if hasattr(self,'eta_v'):
            print('    Volumetric:', self.eta_v,'-')
        
        #Want to return a list
        return [Qnet]
    
    def S_to_SS(self, FlowPath, **kwargs):
        if self.__before_discharge1__:
            FlowPath.A = 0.0
        else:
            FlowPath.A = scroll_geo.Area_d_dd(2*pi-self.beta, self.geo)
        try:
            mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                               FlowPath.State_up,
                                               FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
        
    def DA_D(self, FlowPath, X_d=1.0,**kwargs):
        FlowPath.A=X_d*scroll_geo.Area_s_sa(2*pi-self.beta, self.geo)
        try:
            return flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
        except ZeroDivisionError:
            return 0.0
        
    def step_callback(self,t,h,Itheta):
        """
        Here we test whether the control volumes need to be
        a) Split
        b) Adjusted because you are at the discharge angle
        
        """ 
        #This gets called at every step, or partial step
        self.beta=t
        
        def angle_difference(angle1,angle2):
            # Due to the periodicity of angles, you need to handle the case where the
            # angles wrap around - suppose theta_d is 6.28 and you are at an angles of 0.1 rad
            #, the difference should be around 0.1, not -6.27
            # 
            # This brilliant method is from http://blog.lexique-du-net.com/index.php?post/Calculate-the-real-difference-between-two-angles-keeping-the-sign
            # and the comment of user tk
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
            print('At the discharge angle')
            ########################
            #Reassign chambers
            ########################
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
            
            print('splitting')
            
            if self.__hasLiquid__==False:

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
        
def Expander():
    
    Ref = 'Nitrogen'
    Tinlet = 300
    pinlet = 2000
    poutlet = 500
    
    SE = ScrollExpander()
    SE.set_scroll_geo(100e-6,2.3,0.005,0.005)
    SE.set_disc_geo('2Arc', r2 = 0)
    
    SE.geo.delta_suction_offset = 0.0e-3
    SE.geo.phi_ie_offset = 0.0
    
    SE.omega = 60*2*pi
    SE.Tamb = 298.0
    
    #Temporarily set the bearing dimensions
    SE.mech = core.struct()
    SE.mech.D_upper_bearing = 0.02540
    SE.mech.L_upper_bearing = 0.031
    SE.mech.c_upper_bearing = 10e-6 #??
    SE.mech.D_crank_bearing = 0.0254
    SE.mech.L_crank_bearing = 0.023
    SE.mech.c_crank_bearing = 10e-6 #??
    SE.mech.D_lower_bearing = 0.01895
    SE.mech.L_lower_bearing = 0.023
    SE.mech.c_lower_bearing = 10e-6 #??
    SE.mech.thrust_ID = 0.05
    SE.mech.thrust_friction_coefficient = 0.02 #Tuned for point A1
    SE.mech.L_ratio_bearings = 5
    SE.mech.mu_oil = 0.008
    SE.mech.scroll_plate_thickness = 0.007
    SE.mech.scroll_density = 8100
    SE.mech.scroll_added_mass = 0
    SE.mech.scroll_plate_diameter = 0.092 #[m]
    m, zcm = SE.calculate_scroll_mass()
    SE.mech.orbiting_scroll_mass = m
    SE.mech.scroll_zcm__thrust_surface = zcm
    
    SE.mech.detailed_analysis = False
    
    SE.h_shell = 10
    SE.A_shell = 0.05
    SE.HTC = 0.5
    
    SE.generator = Motor()
    SE.generator.set_eta(1.0)
    SE.generator.suction_fraction = 0.0
    
    inletState = CPState(Ref,{'T':Tinlet,'P':pinlet})
    T2s = SE.guess_outlet_temp(inletState,poutlet)
    outletState = CPState(Ref,{'T':T2s,'P':poutlet})
    
    SE.geo.delta_radial = 0e-6
    SE.geo.delta_flank = 0e-6

    mdot_guess = inletState.rho*SE.Vdisp/SE.Vratio*SE.omega/(2*pi)
    
    SE.add_tube(Tube(key1 = 'inlet.1',
                     key2 = 'inlet.2',
                     L = 0.03,
                     ID = 0.02,
                     mdot = mdot_guess, 
                     State1 = inletState.copy(),
                     fixed = 1,
                     TubeFcn = SE.TubeCode))
    SE.add_tube(Tube(key1 = 'outlet.1',
                     key2 = 'outlet.2',
                     L = 0.03,
                     ID = 0.02,
                     mdot = mdot_guess, 
                     State2 = outletState.copy(),
                     fixed = 2,
                     TubeFcn = SE.TubeCode))
    
    SE.auto_add_CVs(inletState, outletState)
#    SE.auto_add_leakage(flankFunc = SE.FlankLeakage, 
#                        radialFunc = SE.RadialLeakage)

    FP = FlowPath(key1='outlet.1', 
                  key2='da',
                  MdotFcn=IsentropicNozzleWrapper(),
                  )
    FP.A = pi*0.02**2/4
    SE.add_flow(FP)
    
    SE.add_flow(FlowPath(key1 = 'da', 
                         key2 = 'd1',
                         MdotFcn = SE.DA_D,
                         MdotFcn_kwargs = dict(X_d = 0.5)
                         )
                )
    SE.add_flow(FlowPath(key1 = 'da',
                         key2 = 'd2',
                         MdotFcn = SE.DA_D,
                         MdotFcn_kwargs = dict(X_d = 0.5)
                         )
                )
    
    FP = FlowPath(key1='inlet.2', 
                  key2='ss', 
                  MdotFcn=IsentropicNozzleWrapper(),
                  )
    FP.A = pi*0.007**2/4
    SE.add_flow(FP)
    
    FP = FlowPath(key1='inlet.2', 
                  key2='sss', 
                  MdotFcn=IsentropicNozzleWrapper(),
                  )
    FP.A = pi*0.007**2/4
    SE.add_flow(FP)
    
    SE.add_flow(FlowPath(key1='s1',
                         key2='ss',
                         MdotFcn=SE.S_to_SS))
    SE.add_flow(FlowPath(key1='s2',
                         key2='ss',
                         MdotFcn=SE.S_to_SS))
    
    SE.connect_callbacks(step_callback = SE.step_callback,
                         lumps_energy_balance_callback = SE.lump_energy_balance_callback,
                         endcycle_callback= SE.endcycle_callback)
    
#    beta = np.linspace(0,2*pi,1000)
#    V_e11 = np.array([SE.V_e1(b,1)[0] for b in beta])
#    V_e21 = np.array([SE.V_e2(b,1)[0] for b in beta])
#    V_e12 = np.array([SE.V_e1(b,2)[0] for b in beta])
#    V_e22 = np.array([SE.V_e2(b,2)[0] for b in beta])
#    V_d1 = np.array([SE.V_d1(b)[0] for b in beta])
#    V_d2 = np.array([SE.V_d2(b)[0] for b in beta])
#    V_s1 = np.array([SE.V_s1(b)[0] for b in beta])
#    
#    plt.plot(beta,V_e11)
#    plt.plot(beta,V_e12)
#    plt.plot(beta,V_d1)
#    plt.plot(beta,V_s1)
#    plt.gca().axvline(2*pi-SE.theta_d)
#    plt.show()
    #try:
    SE.eps_cycle = 0.003
    SE.precond_solve(key_inlet='inlet.1',
                     key_outlet='outlet.2',
                     solver_method='RK45',
                     UseNR = False, #If True, use Newton-Raphson ND solver to determine the initial state
                     OneCycle = False,
                     plot_every_cycle= False,
                     hmin = 1e-8
                     )
    #except ValueError as E:
    #    print E
    #    pass
    
if __name__=='__main__':
    Expander()
