from __future__ import division
##--- package imports
from PDSim.core.containers import ControlVolume
from PDSim.flow.flow import FlowPath
from PDSim.core.core import PDSimCore
from PDSim.misc._listmath import listm
from PDSim.flow import flow_models
from PDSim.plot.plots import debug_plots
from PDSim.core.bearings import journal_bearing
import scroll_geo
from _scroll import _Scroll

##--- non-package imports
import warnings
from scipy.optimize import fsolve, newton
from CoolProp.CoolProp import UseSinglePhaseLUT, Props
from CoolProp import State
from math import pi,cos
import numpy as np
import copy

class struct(object):
    pass

class Scroll(PDSimCore, _Scroll):
    """
    This is a python class that implements functionality for a scroll compressor
    
    It is inherited from the PDSimCore class
    """
    def __init__(self):
        PDSimCore.__init__(self)
        
        ## Define the geometry structure
        self.geo=scroll_geo.geoVals()
        
        ## Set flags
        self.__Setscroll_geo__=False
        self.__SetDiscGeo__=False
        self.__before_discharge1__=False #Step bridging theta_d
        self.__before_discharge2__=False #Step up to theta_d

    def __getstate__(self):
        """
        A function for preparing class instance for pickling
         
        Combine the dictionaries from the _Scroll base class and the Scroll
        class when pickling
        """
        py_dict = self.__dict__.copy()
        py_dict.update(_Scroll.__cdict__(self))
        return py_dict

    def __setstate__(self, d):
        """
        A function for unpacking class instance for unpickling
        """
        for k,v in d.iteritems():
            setattr(self,k,v)
        
    @property
    def theta_d(self):
        return scroll_geo.theta_d(self.geo)
    
    @property
    def Vdisp(self):
        return -2*pi*self.geo.h*self.geo.rb*self.geo.ro*(3*pi
                                                         -2*self.geo.phi_ie
                                                         +self.geo.phi_i0
                                                         +self.geo.phi_o0)
    
    @property
    def Vratio(self):
        return ((3*pi-2*self.geo.phi_ie+self.geo.phi_i0+self.geo.phi_o0)
                /(-2*self.geo.phi_os-3*pi+self.geo.phi_i0+self.geo.phi_o0))
    
    def V_injection(self, theta, V_tube = None):
        """
        Volume code for injection tube
        
        The tube volume can either be given by the keyword argument V_tube 
        (so you can easily have more than one injection tube), or it can be 
        provided by setting the Scroll class attribute V_inj_tube 
        and NOT providing the V_tube argument
        
        The injection tube volume is assumed to be constant, hence the derivative of volume is zero 
        """
        if V_tube is None:
            return self.V_inj_tube, 0.0
        else:
            return V_tube, 0.0
        
    def V_sa(self, theta, full_output=False):
        """
        Wrapper around the Cython code for sa calcs
        
        Parameters
        ----------
        theta: float
             angle in range [0,2*pi]
        
        Returns
        -------
        
        """
        return scroll_geo.SA(theta,self.geo)[0:2]
        
    def V_s1(self,theta):
        """
        Wrapper around the Cython code for Vs1_calcs
        
        theta: angle in range [0,2*pi]
        """
        return scroll_geo.S1(theta,self.geo)[0:2]
        
    def V_s2(self,theta):
        """
        Wrapper around the Cython code for Vs1_calcs
        
        theta: angle in range [0,2*pi]
        """

        return scroll_geo.S2(theta,self.geo)[0:2]
    
    def V_c1(self,theta,alpha=1,full_output=False):
        """
        Wrapper around the Cython code for C1
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        return scroll_geo.C1(theta,alpha,self.geo)[0:2]
        
    def V_c2(self,theta,alpha=1,full_output=False):
        """
        Wrapper around the Cython code for C2
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        return scroll_geo.C2(theta,alpha,self.geo)[0:2]
        
    def V_d1(self,theta,full_output=False):
        """
        Wrapper around the compiled code for D1
        
        theta: angle in range [0,2*pi]
        """
        
        if self.__before_discharge1__==True and theta<self.theta_d:
                #Get the number of compression chambers in existence
                alpha=scroll_geo.getNc(theta,self.geo)
                #Use the innermost compression chamber 
                return scroll_geo.C1(theta,alpha,self.geo)[0:2]
        else:
            return scroll_geo.D1(theta,self.geo)[0:2]
    
    def V_d2(self,theta,full_output=False):
        """
        Wrapper around the compiled code for D2
        
        theta: angle in range [0,2*pi]
        """

        if self.__before_discharge1__==True and theta<self.theta_d:
                #Get the number of compression chambers in existence
                alpha=scroll_geo.getNc(theta,self.geo)
                #Use the innermost compression chamber 
                return scroll_geo.C2(theta,alpha,self.geo)[0:2]
        else:
            return scroll_geo.D2(theta,self.geo)[0:2]
    
    def V_dd(self,theta,full_output=False):
        """
        Wrapper around the compiled code for DD
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return scroll_geo.DD(theta,self.geo)[0:2],HTangles
        else:
            if self.__before_discharge1__==True and theta<self.theta_d:
                return scroll_geo.DDD(theta,self.geo)[0:2]
            else:
                return scroll_geo.DD(theta,self.geo)[0:2]
        
    def V_ddd(self,theta,alpha=1,full_output=False):
        """
        Wrapper around the compiled code for DDD
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return scroll_geo.DDD(theta,self.geo)[0:2],HTangles
        else:
            return scroll_geo.DDD(theta,self.geo)[0:2]        
        
    def set_scroll_geo(self,Vdisp,Vratio,Thickness,OrbitingRadius,phi_i0=0.0,phi_os=0.3, phi_is = pi):
        """
        Provide the following parameters.  The rest will be calculated by the geometry code
        
        ==============  ===================================================================================
        Vdisp           Displacement in compressor mode [m^3]
        Vratio          Volume ratio (compression chambers at discharge angle / displacement volume) [-]
        Thickness       Thickness of scroll wrap [m]
        OrbitingRadius  Orbiting radius of the orbiting scroll [m]
        ==============  ===================================================================================
        
        Optional parameters are 
        
        phi_i0
        phi_os
        phi_is
        """
        
        ## Determine the geometry by using the imposed parameters for the scroll wraps
        def f(x,phi_i0,phi_os,Vdisp_goal,Vratio_goal,t_goal,ro_goal):
            phi_ie=x[0]
            phi_o0=x[1]
            hs=x[2]
            rb=x[3]
            t=rb*(phi_i0-phi_o0)
            ro=rb*pi-t
            Vdisp=-2*pi*hs*rb*ro*(3*pi-2*phi_ie+phi_i0+phi_o0)
            Vratio=(3*pi-2*phi_ie+phi_i0+phi_o0)/(-2*phi_os-3*pi+phi_i0+phi_o0)

            r1=Vdisp-Vdisp_goal
            r2=Vratio-Vratio_goal
            r3=t-t_goal
            r4=ro-ro_goal
            return [r1,r2,r3,r4]
        
        phi_ie,phi_o0,hs,rb=fsolve(f,[20,1.3,0.03,0.003],args=(phi_i0,phi_os,Vdisp,Vratio,Thickness,OrbitingRadius))
        phi_oe=phi_ie
        self.geo.h=hs
        self.geo.rb=rb
        self.geo.phi_i0=phi_i0
        self.geo.phi_is=phi_is
        self.geo.phi_ie=phi_ie
        self.geo.phi_o0=phi_o0
        self.geo.phi_os=phi_os
        self.geo.phi_oe=phi_oe
        self.geo.ro=rb*pi-Thickness
        self.geo.t=Thickness
        
        #Set the flags to ensure all parameters are fresh
        self.__Setscroll_geo__=True
        self.__SetDiscGeo__=False     
    
    def set_disc_geo(self,Type,r2=0.0):
        """
        Set the discharge geometry for the scrolls
        
        Parameters
        ----------
        Type
            The type of 
        """
        if self.__Setscroll_geo__==False:
            raise ValueError("You must determine scroll wrap geometry by calling Setscroll_geo before setting discharge geometry.")
        
        #Use the compiled code
        scroll_geo.setDiscGeo(self.geo,Type,r2)
        
    def auto_add_CVs(self,inletState,outletState):
        """
        Adds all the control volumes for the scroll compressor.
        
        Parameters
        ----------
        inletState
            A ``State`` instance for the inlet to the scroll set.  Can be approximate
        outletState
            A ``State`` instance for the outlet to the scroll set.  Can be approximate
            
        Notes
        -----
        Uses the indices of 
        
        ============= ===================================================================
        CV            Description
        ============= ===================================================================
        ``sa``        Suction Area
        ``s1``        Suction chamber on side 1
        ``s2``        Suction chamber on side 2
        ``d1``        Discharge chamber on side 1
        ``d2``        Discharge chamber on side 2
        ``dd``        Central discharge chamber
        ``ddd``       Merged discharge chamber
        ``c1.i``      The i-th compression chamber on side 1 (i=1 for outermost chamber)
        ``c2.i``      The i-th compression chamber on side 2 (i=1 for outermost chamber)
        ============= ===================================================================
        """
        
        #Add all the control volumes that are easy.  Suction area and suction chambera
        self.add_CV(ControlVolume(key='sa',initialState=inletState.copy(),
                                VdVFcn=self.V_sa,becomes=['sa','s1','s2']))
        self.add_CV(ControlVolume(key='s1',initialState=inletState.copy(),
                                VdVFcn=self.V_s1,becomes='c1.1'))
        self.add_CV(ControlVolume(key='s2',initialState=inletState.copy(),
                                VdVFcn=self.V_s2,becomes='c2.1'))
        
        #Discharge chambers are also easy.  Assume that you start with 'ddd' chamber merged.
        # No problem if this isn't true.
        self.add_CV(ControlVolume(key='d1',initialState=outletState.copy(),
                                VdVFcn=self.V_d1,exists=False))
        self.add_CV(ControlVolume(key='d2',initialState=outletState.copy(),
                                VdVFcn=self.V_d2,exists=False))
        self.add_CV(ControlVolume(key='dd',initialState=outletState.copy(),
                                VdVFcn=self.V_dd,exists=False))
        self.add_CV(ControlVolume(key='ddd',initialState=outletState.copy(),
                                VdVFcn=self.V_ddd,discharge_becomes='dd'))

        #Add each pair of compression chambers
        nCmax = scroll_geo.nC_Max(self.geo)
        # Must have at least one pair
        assert (nCmax>=1)
        for alpha in range(1,nCmax+1):
            keyc1 = 'c1.'+str(alpha)
            keyc2 = 'c2.'+str(alpha)
            if alpha==1:
                #It is the outermost pair of compression chambers
                initState = State.State(inletState.Fluid,
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
                V1 = self.V_s1(2*pi)[0]
                V2 = self.V_c1(0,alpha)[0]
                #Mass is constant, so rho1*V1 = rho2*V2
                rho2 = rho1 * V1 / V2
                # Now don't know temperature or pressure, but you can assume
                # it is isentropic to find the temperature
                T2 = newton(lambda T: Props('S','T',T,'D',rho2,inletState.Fluid)-s1, T1)
                initState=State.State(inletState.Fluid,dict(T=T2,D=rho2)).copy()
            if alpha<nCmax:
                # Does not change definition at discharge angle
                disc_becomes_c1 = 'c1.'+str(alpha)
                disc_becomes_c2 = 'c2.'+str(alpha)
                # It is not the innermost pair of chambers, becomes another 
                # set of compression chambers at the end of the rotation
                becomes_c1 = 'c1.'+str(alpha+1)
                becomes_c2 = 'c2.'+str(alpha+1)
            else:
                #It is the innermost pair of chambers, becomes discharge chamber
                #at the discharge angle
                disc_becomes_c1 = 'd1'
                disc_becomes_c2 = 'd2'
                becomes_c1 = 'c1.'+str(alpha+1) #Not used - CV dies at disc.
                becomes_c2 = 'c2.'+str(alpha+1) #Not used - CV dies at disc.
                
            self.add_CV(ControlVolume(key=keyc1,
                                      initialState=initState.copy(),
                                      VdVFcn=self.V_c1,
                                      VdVFcn_kwargs={'alpha':alpha},
                                      discharge_becomes=disc_becomes_c1,
                                      becomes=becomes_c1))
            
            self.add_CV(ControlVolume(key=keyc2,
                                      initialState=initState.copy(),
                                      VdVFcn=self.V_c2,
                                      VdVFcn_kwargs={'alpha':alpha},
                                      discharge_becomes=disc_becomes_c2,
                                      becomes=becomes_c2))
    
    def auto_add_leakage(self,flankFunc,radialFunc):
        """
        Add all the leakage terms for the compressor
        
        Parameters
        ----------
        flankFunc : function
            The function to be used for the flank leakage path
        radialFunc : function
            The function to be used for the radial leakage path
        """
        
        #Do the flank leakages
        self.auto_add_flank_leakage(flankFunc)
        #Do the radial leakages
        self.auto_add_radial_leakage(radialFunc)
        
    def auto_add_radial_leakage(self, radialFunc):
        """
        A function to add all the radial leakage terms
        
        Parameters
        ----------
        radialFunc : function
            The function that will be called for each radial leakage
        """
        #Get all the radial leakage pairs
        pairs = scroll_geo.radial_leakage_pairs(self.geo)
        
        #Loop over all the radial leakage pairs possible for the given geometry
        for pair in pairs:
            self.add_flow(FlowPath(key1=pair[0],
                                   key2=pair[1],
                                   MdotFcn=radialFunc
                                   )
                          )
        
    def auto_add_flank_leakage(self, flankFunc):
        """
        A function to add all the flank leakage terms
        
        Parameters
        ----------
        flankFunc : function
            The function that will be called for each flank leakage
        """
        
        # Always a s1-c1 leakage and s2-c2 leakage
        self.add_flow(FlowPath(key1='s1',key2='c1.1',MdotFcn=flankFunc))
        self.add_flow(FlowPath(key1='s2',key2='c2.1',MdotFcn=flankFunc))
        
        # Only add the DDD-S1 and DDD-S2 flow path if there is one set of
        # compression chambers.   
        if scroll_geo.nC_Max(self.geo) == 1:
            self.add_flow(FlowPath(key1='s1',key2='ddd',MdotFcn=self.DDD_to_S))
            self.add_flow(FlowPath(key1='s2',key2='ddd',MdotFcn=self.DDD_to_S))
        
        #Add each pair of compression chambers
        nCmax = scroll_geo.nC_Max(self.geo)
        
        # Must have at least one pair
        assert (nCmax>=1)
        
        
        for alpha in range(1,nCmax+1):
            keyc1 = 'c1.'+str(alpha)
            keyc2 = 'c2.'+str(alpha)
            
            if alpha < nCmax - 1:
                #Leakage between compression chambers along a path
                self.add_flow(FlowPath(key1=keyc1,
                                       key2='c1.'+str(alpha+1),
                                       MdotFcn=flankFunc))
                self.add_flow(FlowPath(key1=keyc2,
                                       key2='c2.'+str(alpha+1),
                                       MdotFcn=flankFunc))
                
            elif alpha==nCmax:
                #Leakage between the discharge region and the innermost chamber
                self.add_flow(FlowPath(key1=keyc1,key2='ddd',MdotFcn=flankFunc))
                self.add_flow(FlowPath(key1=keyc2,key2='ddd',MdotFcn=flankFunc))
    
    def heat_transfer_coefficient(self, key):
        
#        Pr=Pr_mix(Ref,Liq,T_avg,p_avg,xL_avg); //[-]
#        Re=4.0*mdot/2.0/(PI*mu_mix(Ref,Liq,T_avg,p_avg,xL_avg)*Dh); //[-]
#        hc=0.023*k_mix(Ref,Liq,T_avg,p_avg,xL_avg)/Dh*pow(Re,0.8)*pow(Pr,0.4); //[kW/m^2-K]
#        // Jang and Jeong correction for spiral geometry
#        f=scroll->States.omega/(2*PI);
#        Amax=scroll->geo.ro;
#        Ubar=scroll->massFlow.mdot_tot/(4*scroll->geo.ro*scroll->geo.hs*rho);
#        St=f*Amax/Ubar;
#        hc*=1.0+8.48*(1-exp(-5.35*St));
#        // Tagri and Jayaraman correction for transverse oscillation
#        r_c=scroll->geo.rb*(0.5*phi_1_i+0.5*phi_2_i-scroll->geo.phi.phi_fi0);
#        hc*=1.0+1.77*Dh/r_c;
        return 1.0

    def wrap_heat_transfer(self, **kwargs):
        """
        This function evaluates the anti-derivative of 
        the differential of wall heat transfer, and returns the amount of scroll-
        wall heat transfer in kW
        
        Parameters
        ----------
        hc : float
            Heat transfer coefficient [kW/m2/K]
        hs : float
            Scroll wrap height [m]
        rb : float
            Base circle radius [m]
        phi1 : float
            Larger involute angle [rad]
        phi2 : float
            Smaller involute angle [rad]
        phi0 : float
            Initial involute angle [rad]
        T_scroll : float
            Lump temperature of the scroll wrap [K]
        T_CV : float
            Temperature of the gas in the CV [K]
        dT_dphi : float
            Derivative of the temperature along the scroll wrap [K/rad]
        phim : float
            Mean involute angle of wrap used for heat transfer [rad]
        
        Notes
        -----
        ``phi1`` and ``phi2`` are defined such that ``phi1`` is always the
        larger involute angle in value
        """
        #Use the compiled version from the cython code
        return _Scroll.involute_heat_transfer(self,**kwargs)
    
    def heat_transfer_callback(self, theta, **kwargs):
        """
        The scroll simulation heat transfer callback for HT to the fluid in the 
        chambers
        
        ``heat_transfer_callback`` for ``PDSimCore.derivs`` must be of the 
        form::
        
            heat_transfer_callback(theta, **kwargs)
            
        but we need to get the inlet and outlet states to get the linear 
        temperature profile in the scroll wrap. Thus we wrap the callback 
        we would like to call in this function that allows us to determine
        the inlet and outlet state at run-time.
        """
        State_inlet = self.Tubes.Nodes[self.key_inlet]
        State_outlet = self.Tubes.Nodes[self.key_outlet]
        return self._heat_transfer_callback(theta, State_inlet, State_outlet, **kwargs)
    
    def _heat_transfer_callback(self, theta, State_inlet, State_outlet, HTC_tune = 1.0, **kwargs):
        """
        A private function to actually do the heat transfer analysis
        """
        # dT_dphi is generally negative because as you move to the 
        # outside of the scroll (larger phi), the temperature goes down because
        # you are moving towards low pressure and low temperature

        Tsuction = State_inlet.T
        Tdischarge = State_outlet.T
        dT_dphi = (Tsuction - Tdischarge) / (self.geo.phi_ie - self.geo.phi_os)
        phim = 0.5*self.geo.phi_ie + 0.5*self.geo.phi_os
        
        def calcHT(key):
            try:
                ## If HT is turned off, quit
                if HTC_tune <= 0.0:
                    return 0.0
                    
                #Get the bounding angles for the control volume
                HTangles = scroll_geo.HT_angles(theta, self.geo, key)
                
                phi_1_i=HTangles['1_i']
                phi_2_i=HTangles['2_i']
                phi_1_o=HTangles['1_o']
                phi_2_o=HTangles['2_o']
    
                #print 'calculate HTC'
                #TODO: calculate HTC
                hc = self.HTC #[kW/m2/K]
                
                # The heat transfer rate of the inner involute on 
                # the outer wrap of the chamber
                Q_outer_wrap = self.wrap_heat_transfer(hc = hc, 
                                                   hs = self.geo.h, 
                                                   rb = self.geo.rb, 
                                                   phi1 = phi_1_i, 
                                                   phi2 = phi_2_i, 
                                                   phi0 = self.geo.phi_i0, 
                                                   T_scroll = self.Tlumps[0],
                                                   T_CV = self.CVs[key].State.T, 
                                                   dT_dphi = dT_dphi, 
                                                   phim = phim)
                
                # The heat transfer rate of the outer involute on 
                # the inner wrap of the chamber
                Q_inner_wrap = self.wrap_heat_transfer(hc = hc, 
                                                   hs = self.geo.h, 
                                                   rb = self.geo.rb, 
                                                   phi1 = phi_1_o, 
                                                   phi2 = phi_2_o, 
                                                   phi0 = self.geo.phi_o0,
                                                   T_scroll = self.Tlumps[0],
                                                   T_CV = self.CVs[key].State.T, 
                                                   dT_dphi = dT_dphi, 
                                                   phim = phim)
                
                return HTC_tune *(Q_outer_wrap + Q_inner_wrap)
            
            except KeyError:
                if key == 'sa':
                    # sa is treated as having no heat transfer
                    return 0.0
                elif key == 'ddd':
                    # ddd is a combination of the heat transfer in the d1, d2, and
                    # dd chambers
                    return calcHT('d1') + calcHT('d2') + calcHT('dd')
                elif key == 'dd':
                    # dd chamber is treated as having no heat transfer
                    return 0.0
                elif key.startswith('inj'):
                    # injection chambers are treated as having no heat transfer
                    return 0.0
                else:
                    raise KeyError('CV '+key+' not matched')
            
        Q=map(calcHT, self.CVs.exists_keys)
        return listm(Q)
        
    def step_callback(self,t,h,Itheta,**kwargs):
        """
        Here we test whether the control volumes need to be
        a) Merged
        b) Adjusted because you are at the discharge angle
        
        """ 
        #This gets called at every step, or partial step
        self.theta=t
        
        def IsAtMerge(eps = 0.005, eps_d1_higher=0.005,eps_dd_higher=0.00001):
            pressures = [self.CVs['d1'].State.p,
                         self.CVs['d2'].State.p,
                         self.CVs['dd'].State.p]
            p_max = max(pressures)
            p_min = min(pressures)
            if abs(p_min/p_max-1)<eps_dd_higher:
                return True
            elif self.CVs['d1'].State.p>self.CVs['dd'].State.p and abs(p_min/p_max-1)<eps_d1_higher:
                print 'Merged with d1 higher'
                return True
            elif self.CVs['d1'].State.p<self.CVs['dd'].State.p and abs(p_min/p_max-1)<eps_dd_higher:
                print 'Merged with dd higher'
                return True
            else:
                return False
            
        disable=False
        
        if t<self.theta_d<t+h and self.__before_discharge2__==False:
            #Take a step almost up to the discharge angle
            disable=True
            h=self.theta_d-t-1e-10
            self.__before_discharge2__=True
        elif self.__before_discharge2__==True:
            #At the discharge angle
            print 'At the discharge angle'
            ########################
            #Reassign chambers
            ########################
            #Find chambers with a discharge_becomes flag
            for key in self.CVs.exists_keys:
                if self.CVs[key].discharge_becomes in self.CVs.keys():
                    #Set the state of the "new" chamber to be the old chamber
                    oldCV=self.CVs[key]
                    if oldCV.exists==True:
                        newCV=self.CVs[oldCV.discharge_becomes]
                        newCV.State.update({'T':oldCV.State.T,'D':oldCV.State.rho})
                        oldCV.exists=False
                        newCV.exists=True
                    else:
                        raise AttributeError("old CV doesn't exist")
            
            self.__before_discharge2__=False
            self.__before_discharge1__=True
            
            self.CVs.rebuild_exists()
            #Re-calculate the CV volumes
            V,dV = self.CVs.volumes(t)
            #Update the matrices using the new CV definitions
            self.T[self.CVs.exists_indices,Itheta]=self.CVs.T
            self.p[self.CVs.exists_indices,Itheta]=self.CVs.p
            self.m[self.CVs.exists_indices,Itheta]=listm(self.CVs.rho)*V
            self.rho[self.CVs.exists_indices,Itheta]=listm(self.CVs.rho)
            
            # Adaptive makes steps of h/4 3h/8 12h/13 and h/2 and h
            # Make sure step does not hit any *right* at theta_d
            # That is why it is 2.2e-8 rather than 2.0e-8
            h=2.2e-10
            disable=True
       
        elif self.CVs['d1'].exists and IsAtMerge():
            
            #Build the volume vector using the old set of control volumes (pre-merge)
            V,dV=self.CVs.volumes(t)
            
            if self.__hasLiquid__==False:

                #Density
                rhod1=self.CVs['d1'].State.rho
                rhod2=self.CVs['d2'].State.rho
                rhodd=self.CVs['dd'].State.rho
                #Density
                pd1=self.CVs['d1'].State.p
                pd2=self.CVs['d2'].State.p
                pdd=self.CVs['dd'].State.p
                #Internal energy
                ud1=self.CVs['d1'].State.u
                ud2=self.CVs['d2'].State.u
                udd=self.CVs['dd'].State.u
                #Internal energy
                Td1=self.CVs['d1'].State.T
                Td2=self.CVs['d2'].State.T
                Tdd=self.CVs['dd'].State.T
                #Volumes
                Vdict=dict(zip(self.CVs.exists_keys,V))
                Vd1=Vdict['d1']
                Vd2=Vdict['d2']
                Vdd=Vdict['dd']
                
                Vddd=Vd1+Vd2+Vdd
                m=rhod1*Vd1+rhod2*Vd2+rhodd*Vdd
                U_before=ud1*rhod1*Vd1+ud2*rhod2*Vd2+udd*rhodd*Vdd
                rhoddd=m/Vddd
                #guess the mixed temperature as a volume-weighted average
                T=(Td1*Vd1+Td2*Vd2+Tdd*Vdd)/Vddd
                p=(pd1*Vd1+pd2*Vd2+pdd*Vdd)/Vddd
                #Must conserve mass and internal energy (instantaneous mixing process)
                Fluid = self.CVs['ddd'].State.Fluid
                T_u = newton(lambda x: Props('U','T',x,'D',rhoddd,Fluid)-U_before/m,T)
                
                self.CVs['ddd'].State.update({'T':T_u,'D':rhoddd})
                U_after=self.CVs['ddd'].State.u*self.CVs['ddd'].State.rho*Vddd
                
                DeltaU=m*(U_before-U_after)
                if abs(DeltaU)>1e-5:
                    raise ValueError('Internal energy not sufficiently conserved in merging process')
                
                self.CVs['d1'].exists=False
                self.CVs['d2'].exists=False
                self.CVs['dd'].exists=False
                self.CVs['ddd'].exists=True
                
                self.CVs.rebuild_exists()
                
                #Re-calculate the CV
                V,dV=self.CVs.volumes(t)
                self.T[self.CVs.exists_indices,Itheta]=self.CVs.T
                self.p[self.CVs.exists_indices,Itheta]=self.CVs.p
                self.m[self.CVs.exists_indices,Itheta]=listm(self.CVs.rho)*V
                self.rho[self.CVs.exists_indices,Itheta]=listm(self.CVs.rho)
                
            else:
                raise NotImplementedError('no flooding yet')
            disable=True 
              
        elif t>self.theta_d:
            self.__before_discharge1__=False
            disable=False
            
        return disable,h
        
    def crank_bearing(self):
        
        JB = journal_bearing(r_b = self.D_crank_bearing/2,
                             L = self.L_crank_bearing,
                             omega = self.omega,
                             W = self.forces.corr_Fm*1000,
                             c = self.c_crank_bearing,
                             eta_0 = self.mu_oil
                             )
    
        #print 'Crank pin journal loss is',JB['Wdot_loss'],'W'
        return JB['Wdot_loss']/1000.0
        
    def upper_bearing(self):
        """
        Moment balance around the upper bearing gives the force for
        the lower bearing.  Torques need to balance around the upper bearing
        """
        
        JB = journal_bearing(r_b = self.D_upper_bearing/2,
                             L = self.L_upper_bearing,
                             omega = self.omega,
                             W = self.forces.corr_Fm*1000*(1+1/self.L_ratio_bearings),
                             c = self.c_upper_bearing,
                             eta_0 = self.mu_oil
                             )
        #print 'Upper bearing journal loss is',JB['Wdot_loss'],'W'
        return JB['Wdot_loss']/1000.0
    
    def lower_bearing(self):
        """
        Moment balance around the upper bearing gives the force for
        the lower bearing.  Torques need to balance around the upper bearing
        """
        
        JB = journal_bearing(r_b = self.D_lower_bearing/2,
                             L = self.L_lower_bearing,
                             omega = self.omega,
                             W = self.forces.corr_Fm*1000*(1/self.L_ratio_bearings),
                             c = self.c_lower_bearing,
                             eta_0 = self.mu_oil
                             )
        
        #print 'Lower bearing journal loss is',JB['Wdot_loss'],'W'
        return JB['Wdot_loss']/1000.0
    
    def thrust_bearing(self):
        """
        The thrust bearing analysis
        """
        from PDSim.core.bearings import thrust_bearing
        V = self.geo.ro*self.omega
        #Use the corrected force to account for the decrease in back area due to the bearing
        N = self.forces.corr_Fz*1000 #[N]
        TB = thrust_bearing(mu = self.thrust_friction_coefficient,
                            V = V,
                            N = N)
        return TB['Wdot_loss']/1000.0
    
    def mechanical_losses(self, shell_pressure = 'low'):
        """
        Calculate the mechanical losses in the bearings
        
        Parameters
        ----------
            shell_pressure : string, 'low' or 'high'

        """
        
        #inlet pressure [kPa]
        inlet_pressure = self.Tubes.Nodes[self.key_inlet].p
        outlet_pressure = self.Tubes.Nodes[self.key_outlet].p
        
        # Get the shell pressure based on either the inlet or outlet pressure
        # based on whether it is a low-pressure or high-pressure shell
        if shell_pressure == 'low':
            back_pressure = min((inlet_pressure, outlet_pressure))
        elif shell_pressure == 'high':
            back_pressure = max((inlet_pressure, outlet_pressure))
        elif shell_pressure == 'mid':
            back_pressure = (inlet_pressure + outlet_pressure)/2
        else:
            raise KeyError("keyword argument shell_pressure must be one of 'low', 'mid' or 'high'")
        
        #Calculate the force terms: force profiles, mean values, etc. 
        self.calculate_force_terms(orbiting_back_pressure = back_pressure)
        
        if not hasattr(self,'losses'):
            self.losses = struct()
            
        if not hasattr(self,'journal_tune_factor'):
            self.journal_tune_factor = 1.0
            
        #Conduct the calculations for the bearings
        self.losses.crank_bearing = self.crank_bearing()*self.journal_tune_factor
        self.losses.upper_bearing = self.upper_bearing()*self.journal_tune_factor
        self.losses.lower_bearing = self.lower_bearing()*self.journal_tune_factor
        self.losses.thrust_bearing = self.thrust_bearing()
        
        self.losses.bearings  = (self.losses.crank_bearing 
                                 + self.losses.upper_bearing 
                                 + self.losses.lower_bearing
                                 + self.losses.thrust_bearing)
        
        print 'mechanical losses: ', self.losses.bearings
        return self.losses.bearings #[kW]
    
    def post_solve(self):
        """
        """
        self.mechanical_losses('low')
        print 'mean_Q', self.HTProcessed.mean_Q
        PDSimCore.post_solve(self)
        
    def ambient_heat_transfer(self, Tshell):
        """
        The amount of heat transfer from the compressor to the ambient
        """
        return self.h_shell*self.A_shell*(Tshell-self.Tamb)
    
    def initial_motor_losses(self, eta_a = 0.8):
        """
        Assume a 70% adiabatic efficiency to estimate the motor power and 
        motor losses
        """
        
        for Tube in self.Tubes:
            if self.key_inlet in [Tube.key1, Tube.key2]:
                mdot = Tube.mdot
                
        inletState = self.Tubes.Nodes[self.key_inlet]
        outletState = self.Tubes.Nodes[self.key_outlet]
        s1 = inletState.s
        h1 = inletState.h
        h2s = Props('H', 'S', s1, 'P', outletState.p, inletState.Fluid)
        
        if outletState.p > inletState.p:
            #Compressor Mode
            h2 = h1 + (h2s-h1)/eta_a
        else:
            #Expander Mode
            h2 = h1 + (h2s-h1)*eta_a
        
        # A guess for the compressor mechanical power based on 70% efficiency [kW]
        Wdot = abs(mdot*(h2-h1))
        
        if self.motor.type == 'const_eta_motor':
            eta = self.motor.eta_motor
        else:
            #The efficiency and speed [-,rad/s] from the mechanical power output
            eta, self.omega = self.motor.invert_map(Wdot)
        
        #Motor losses [kW]
        self.motor.losses = Wdot*(1.0/eta-1)
        
    def suction_heating(self):
        if hasattr(self,'motor'):
            # If some fraction of heat from motor losses is going to get added
            # to suction flow
            if 0.0 <= self.motor.suction_fraction <= 1.0:
                for Tube in self.Tubes:
                    # Find the tube that has one of the keys starting with 'inlet'
                    if Tube.key1.startswith('inlet') or Tube.key2.startswith('inlet'):
                        #Add some fraction of the motor losses to the inlet gas 
                        Tube.Q_add = self.motor.losses * self.motor.suction_fraction
                    else:
                        Tube.Q_add = 0.0
                        
    def pre_run(self):
        """
        Intercepts the call to pre_run and does some scroll processing, then 
        calls the base class function
        """
        
        #Get an initial guess before running at all for the motor losses.
        self.initial_motor_losses()
        
        #Run the suction heating code
        self.suction_heating()
        
        #Call the base class function        
        PDSimCore.pre_run(self)
        
        
    def lump_energy_balance_callback(self):
        """
        .. math ::
            
            \\eta _{motor}} = \\frac{{{{\dot W}_{shaft}}}}{{{{\dot W}_{shaft}} + {{\dot W}_{motor}}}}
            
        .. math ::
            
            {\eta _{motor}}\left( {{{\dot W}_{shaft}} + {{\dot W}_{motor}}} \\right) = {{\dot W}_{shaft}}
            
        .. math::
        
            {{\dot W}_{motor}} = \frac{{{{\dot W}_{shaft}}}}{{{\eta _{motor}}}} - {{\dot W}_{shaft}}
        """
        
        #For the single lump
        # HT terms are positive if heat transfer is TO the lump
        Qnet=0.0
        for Tube in self.Tubes:
            Qnet-=Tube.Q
        
        self.Qamb = self.ambient_heat_transfer(self.Tlumps[0])
        
        # Heat transfer with the ambient; Qamb is positive if heat is being removed, thus flip the sign
        Qnet -= self.Qamb
        
        Qnet += self.mechanical_losses('low') 
        # Heat transfer with the gas in the working chambers.  mean_Q is positive
        # if heat is transfered to the gas in the working chamber, so flip the 
        # sign for the lump
        Qnet -= self.HTProcessed.mean_Q
        
        print 'mean_Q', self.HTProcessed.mean_Q
        
        #Shaft power from forces on the orbiting scroll from the gas in the pockets [kW]
        self.Wdot_forces = self.omega*self.forces.mean_tau
        
        self.Wdot_mechanical = self.Wdot_pv + self.losses.bearings
        
        #The actual torque required to do the compression [N-m]
        self.tau_mechanical = self.Wdot_mechanical / self.omega * 1000
        
        # 2 Options for the motor losses:
        # a) Constant efficiency
        # b) Based on the torque-speed-efficiency motor
        
        if self.motor.type == 'const_eta_motor':
            self.eta_motor = self.motor.eta_motor
        elif self.motor.type == 'motor_map':
            # Use the motor map to calculate the slip rotational speed [rad/s]
            # and the motor efficiency as a function of the torque [N-m]
            eta, omega = self.motor.apply_map(self.tau_mechanical)
            self.eta_motor = eta
            self.omega = omega
        else:
            raise AttributeError
        
        print 'self.forces.mean_Fm', self.forces.mean_Fm
        print 'self.forces.inertial', self.forces.inertial
        print 'self.Qamb', self.Qamb
        print 'self.Wdot_forces', self.Wdot_forces
        print 'self.Wdot_pv', self.Wdot_pv 
        print 'self.losses.bearings', self.losses.bearings
        print 'self.Wdot_mechanical', self.Wdot_mechanical
        
        #Motor losses [kW]
        self.motor.losses = self.Wdot_mechanical*(1/self.eta_motor-1)
        
        #Electrical Power
        self.Wdot_electrical = self.Wdot_mechanical + self.motor.losses
        
        #Overall isentropic efficiency
        self.eta_oi = self.Wdot_i/self.Wdot_electrical
        
#        #Set the heat input to the suction line
#        self.suction_heating()
        
        print 'At this iteration'
        print '    Electrical power:', self.Wdot_electrical
        print '    Mass flow rate:', self.mdot
        
        #Want to return a list
        return [Qnet]
    
    def TubeCode(self,Tube,**kwargs):
        Tube.Q = flow_models.IsothermalWallTube(Tube.mdot,
                                                Tube.State1,
                                                Tube.State2,
                                                Tube.fixed,
                                                Tube.L,
                                                Tube.ID,
                                                T_wall=self.Tlumps[0],
                                                Q_add = Tube.Q_add,
                                                alpha = Tube.alpha
                                                )
        
    
    def DDD_to_S(self,FlowPath,flankFunc = None,**kwargs):
        if  flankFunc is None:
            flankFunc = self.FlankLeakage
        # If there are any compression chambers, don't evaluate this flow
        # since the compression chambers "get in the way" of flow directly from 
        # ddd to s1 and s2
        if scroll_geo.getNc(self.theta,self.geo) > 0:
            return 0.0
        else:
            return flankFunc(FlowPath)
            
    def D_to_DD(self,FlowPath,**kwargs):
        if self.__before_discharge1__:
            FlowPath.A = 0.0
        else:
            FlowPath.A=scroll_geo.Area_d_dd(self.theta,self.geo)
        try:
            return flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
        except ZeroDivisionError:
            return 0.0
        
#    def Inlet_sa(self,*args,**kwargs):
#        mdot = _Scroll.Inlet_sa(self,*args,**kwargs)
#        return mdot
#        
#    def Discharge(self,*args,**kwargs):
#        return _Scroll.Discharge(self,*args,**kwargs)
#        
#    def SA_S(self,*args,**kwargs):
#        return _Scroll.SA_S(self,*args,**kwargs)
#        
#    def FlankLeakage(self,*args,**kwargs):
#        return _Scroll.FlankLeakage(self,*args,**kwargs)
     
    def SA_S1(self, FlowPath, X_d=1.0,**kwargs):
        """
        A wrapper for the flow between the suction area and the S1 chamber
        
        Notes
        -----
        If geo.phi_ie_offset is greater than 0, the offset geometry will be 
        used to calculate the flow area.  Otherwise the conventional analysis 
        will be used.
        """
        if self.geo.phi_ie_offset > 0:
            FlowPath.A = X_d*scroll_geo.Area_s_s1_offset(self.theta, self.geo)
        else:
            FlowPath.A = X_d*scroll_geo.Area_s_sa(self.theta, self.geo)
             
        try:
            mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0   
        
    def SA_S2(self, *args, **kwargs):
        """
        A thin wrapper to the default suction area-suction flow
        """
        return self.SA_S(*args,**kwargs)
        
    def SA_S(self, FlowPath, X_d=1.0,**kwargs):
        
        FlowPath.A=X_d*scroll_geo.Area_s_sa(self.theta, self.geo)
        try:
            mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
    
    def RadialLeakage(self,FlowPath,**kwargs):
        """
        Calculate the radial leakge flow rate
        """
        return _Scroll.RadialLeakage(self, FlowPath, kwargs)
    
        #Calculate the area
        #Arc length of the upstream part of the flow path
        try:
            FlowPath.A= scroll_geo.radial_leakage_area(self.theta,
                                                       self.geo,
                                                       FlowPath.key1,
                                                       FlowPath.key2)
        except KeyError:
            print FlowPath.key1,FlowPath.key2,'caused a KeyError'
            return 0.0
            
        try:
            return flow_models.FrictionCorrectedIsentropicNozzle(
                                 FlowPath.A,
                                 FlowPath.State_up,
                                 FlowPath.State_down,
                                 self.geo.delta_radial,
                                 Type = 'radial',
                                 t = self.geo.t
                                 )
        except ZeroDivisionError:
            return 0.0
            
    def FlankLeakage(self,FlowPath,**kwargs):
        """
        Calculate the flank leakge flow rate
        """
        #Calculate the area
        FlowPath.A=self.geo.h*self.geo.delta_flank
        try:
            return flow_models.FrictionCorrectedIsentropicNozzle(
                                 FlowPath.A,
                                 FlowPath.State_up,
                                 FlowPath.State_down,
                                 self.geo.delta_flank,
                                 Type = 'flank',
                                 ro = self.geo.ro
                                 )
        except ZeroDivisionError:
            return 0.0
        
    def _get_injection_CVkey(self,phi,theta,inner_outer):
        """
        Find the CV that is in contact with the given injection port location
        
        Parameters
        ----------
        phi : float
            Involute angle of the injection port location
        theta : float
            Crank angle in radians in the range [:math:`0,2\pi`]
        inner_outer : string ['i','o']
            'i' : involute angle corresponds to outer surface of fixed scroll
            'o' : involute angle corresponds to inner surface of orb. scroll 
            
        Notes
        -----
        Typically 'i' will require a positive offset in involute angle of 
        :math:`\pi` radians
        """
        if inner_outer == 'i':
            phi_0 = self.geo.phi_i0
            phi_s = self.geo.phi_is
            phi_e = self.geo.phi_ie
        elif inner_outer == 'o':
            phi_0 = self.geo.phi_o0
            phi_s = self.geo.phi_os
            phi_e = self.geo.phi_oe-pi # The actual part of the wrap that could 
                                       # have an injection port 
        
        Nc = scroll_geo.getNc(theta, self.geo)    
        #Start at the outside of the given scroll wrap
        # x1 where x is s,d,c has the inner involute of the fixed scroll as 
        # its outer surface
        if phi_e > phi > phi_e-theta:     
            #It is a suction chamber    
            return 's1' if inner_outer == 'i' else 's2'
            
        elif phi_e-theta > phi > phi_e-theta-2*pi*Nc:
            #It is one of the compression chambers, figure out which one
            for I in range(Nc+1):
                if phi_e - theta - 2*pi*(I-1) > phi > phi_e - theta - 2*pi*I:
                    i_str = '.'+str(I)
                    break
            return 'c1'+i_str if inner_outer == 'i' else 'c2'+i_str
        
        else:
            return 'd1' if inner_outer == 'i' else 'd2'
        
    def Injection_to_Comp(self,FlowPath,phi,inner_outer,check_valve = False, **kwargs):
        """
        Function to calculate flow rate between injection line and chamber
        
        Parameters
        ----------
        FlowPath : FlowPath instance
        phi : involute angle where the port is located
        inner_outer : string ['i','o']
            'i' : involute angle corresponds to outer surface of fixed scroll
            'o' : involute angle corresponds to inner surface of orb. scroll 
        check_valve : boolean
            If ``True``, there is an idealized check valve and flow can only go 
            from chambers with key names that start with `injCV` to other chambers.
            If ``False``, flow can go either direction
            
        .. plot::
            import matplotlib.pyplot as plt
            from PDSim.scroll.plots import plotScrollSet
            import numpy as np
            from math import pi
            from PDSim.scroll import scroll_geo
            
            theta = pi/3
            #Plot the scroll wraps
            plotScrollSet(theta, ScrollComp.geo)
            ax = plt.gca()
            
            #Plot the injection ports (symmetric)
            phi = ScrollComp.geo.phi_oe-pi-2*pi
            #Involute angle along the outer involute of the scroll wrap
            x,y = scroll_geo.coords_inv(phi,ScrollComp.geo,theta,'fo')
            nx,ny = scroll_geo.coords_norm(phi,ScrollComp.geo,theta,'fo')
            rport = 0.002
            xc,yc = x-nx*rport,y-ny*rport
            ax.plot(xc, yc, '.')
            t = np.linspace(0,2*pi,100)
            ax.plot(xc + rport*np.cos(t),yc+rport*np.sin(t),'k')
            
            #Plot the injection ports (symmetric)
            phi = ScrollComp.geo.phi_oe-pi-2*pi+pi
            #Involute angle along the outer involute of the scroll wrap
            x,y = scroll_geo.coords_inv(phi,ScrollComp.geo,theta,'fi')
            nx,ny = scroll_geo.coords_norm(phi,ScrollComp.geo,theta,'fi')
            rport = 0.002
            xc,yc = x-nx*rport,y-ny*rport
            ax.plot(xc, yc, '.')
            t = np.linspace(0,2*pi,100)
            ax.plot(xc + rport*np.cos(t),yc+rport*np.sin(t),'k')
            
            plt.show()
        
        """
        #1. Figure out what CV is connected to the port
        partner_key = self._get_injection_CVkey(phi, self.theta, inner_outer)

        #2. Based on what CV is connected to the port, maybe quit
        if partner_key in ['d1', 'd2'] and 'ddd' in [FlowPath.key_up, 
                                                     FlowPath.key_down]:
            # Other chamber based on geometry is d1 or d2 but they are not 
            # defined due to the angle but ddd is, and as a result, use 
            # ddd
            #
            # Don't do anything, just let it go to the next section even though
            # 'd1' or 'd2 is not key_up or key_down 
            pass
        
        elif partner_key not in [FlowPath.key_up, FlowPath.key_down]:
            return 0.0
        # If the pressure in the injection line is below the other chamber and 
        # you are using a theoretical check valve with instantaneous closing, 
        # then you need to 
        elif check_valve and FlowPath.key_down.startswith('injCV'):
            return 0.0
        #3. Find the distance of the scroll from the point on the involute
        #   where the port is tangent
        FlowPath.A = 3*pi*(0.001)**2/4.0
        try:
            mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
        
        
    def calculate_force_terms(self,
                              orbiting_back_pressure=None):
        """
        Calculate the force profiles, mean forces, moments, etc.
        
        Parameters
        ----------
        orbiting_back_pressure : float, or class instance
            If a class instance, must provide a function __call__ that takes as its first input the Scroll class
        
        """
        
        self.forces = struct()
        
        #Get the slice of indices that are in use.  At the end of the simulation
        #execution this will be the full range of the indices, but when used
        # at intermediate iterations it will be a subset of the indices
        _slice = range(self.Itheta)
        
        ####################################################
        #############  Normal force components #############
        ####################################################
        #The force of the in each chamber pushes the orbiting scroll away
        self.forces.Fz = self.p*self.V/self.geo.h
        
        if isinstance(orbiting_back_pressure, float):
            #The back gas pressure on the orbiting scroll pushes the scroll back down
            self.forces.Fz -= 0.4*orbiting_back_pressure*self.V/self.geo.h
            pass
        else:
            raise NotImplementedError('calculate_force_terms must get a float back pressure for now')
        
        # The Isa chamber doesn't contribute to the total force upwards
        # so zero out its elements in the force matrix
        Isa = self.CVs.index('sa')
        self.forces.Fz[Isa,:] = 0.0
        
        #Remove all the NAN placeholders and replace them with zero values
        self.forces.Fz[np.isnan(self.forces.Fz)]=0
        #Sum the terms
        self.forces.summed_Fz = np.sum(self.forces.Fz, axis = 0) #kN    
        #Calculate the mean axial force
        self.forces.mean_Fz = np.trapz(self.forces.summed_Fz[_slice], self.t[_slice])/(2*pi)
        
        # The corrected axial load accounts for the fact that the area of the 
        # thrust bearing does not contribute to compensating for the load
        
        #: Corrected axial load [kN-m]
        self.forces.corr_Fz = self.forces.mean_Fz# + pi*(self.thrust_OD**2-self.thrust_ID**2)/4 * orbiting_back_pressure
        
        
        
        ####################################################
        #############  "Radial" force components ###########
        ####################################################
        self.forces.Fx = np.zeros((self.CVs.N,len(self.t)))
        self.forces.Fy = np.zeros_like(self.forces.Fx)
        self.forces.Fr = np.zeros_like(self.forces.Fx)
        self.forces.fxp = np.zeros_like(self.forces.Fx)
        self.forces.fyp = np.zeros_like(self.forces.Fx)
        # A map of CVkey to function to be called to get force components
        # All functions in this map use the same call signature and are "boring"
        # Each function returns a dictionary of terms
        func_map = dict(sa = scroll_geo.SA_forces,
                        s1 = scroll_geo.S1_forces,
                        s2 = scroll_geo.S2_forces,
                        d1 = scroll_geo.D1_forces,
                        d2 = scroll_geo.D2_forces,
                        dd = scroll_geo.DD_forces,
                        ddd = scroll_geo.DDD_forces
                        )
        for CVkey in self.CVs.keys():
            if CVkey in func_map:
                #Calculate the force components for each crank angle
                #Early bind the function
                func = func_map[CVkey]
                # Calculate the geometric parts for each chamber
                # They are divided by the pressure in the chamber
                geo_components = [func(theta, self.geo) for theta in self.t[_slice]]
            elif CVkey.startswith('c1'):
                #Early bind the function
                func = scroll_geo.C1_forces
                #Get the key for the CV
                alpha = int(CVkey.split('.')[1])
                # Calculate the geometric parts for each chamber
                # They are divided by the pressure in the chamber
                geo_components = [func(theta,alpha,self.geo) for theta in self.t[_slice]]
            elif CVkey.startswith('c2'):
                #Early bind the function
                func = scroll_geo.C2_forces
                #Get the key for the CV
                alpha = int(CVkey.split('.')[1])
                # Calculate the geometric parts for each chamber
                # They are divided by the pressure in the chamber
                geo_components = [func(theta,alpha,self.geo) for theta in self.t[_slice]]
            else:
                geo_components = []
                
            if geo_components:
                I = self.CVs.index(CVkey)
                p = self.p[I,_slice]
                self.forces.fxp[I,_slice] = [comp['fx_p'] for comp in geo_components]
                self.forces.fyp[I,_slice] = [comp['fy_p'] for comp in geo_components]
                self.forces.Fx[I,_slice] = [comp['fx_p'] for comp in geo_components]*p
                self.forces.Fy[I,_slice] = [comp['fy_p'] for comp in geo_components]*p
        
        #Remove all the NAN placeholders
        self.forces.Fx[np.isnan(self.forces.Fx)]=0
        self.forces.Fy[np.isnan(self.forces.Fy)]=0
        #Sum the terms at each crank angle
        self.forces.summed_Fx = np.sum(self.forces.Fx,axis = 0) #kN
        self.forces.summed_Fy = np.sum(self.forces.Fy,axis = 0) #kN
        
        #Position of the pin as a function of crank angle
        self.forces.xpin = self.geo.ro*np.cos(self.geo.phi_ie-pi/2-self.t)
        self.forces.ypin = self.geo.ro*np.sin(self.geo.phi_ie-pi/2-self.t)
        
        #Calculate the radial force on the crank pin at each crank angle
        
        #The radial component magnitude is just the projection of the force onto a vector going from origin to center of orbiting coordinate system 
        self.forces.Fr = (self.forces.xpin*self.forces.summed_Fx+self.forces.ypin*self.forces.summed_Fy)/self.geo.ro
        
        #The tangent component magnitude is just the projection of the force onto a vector going from origin to center of orbiting coordinate system
        
        #Components of the unit vector in the direction of rotation
        x_dot =  np.sin(self.geo.phi_ie-pi/2-self.t)
        y_dot =  -np.cos(self.geo.phi_ie-pi/2-self.t)
        self.forces.Ft = (x_dot*self.forces.summed_Fx+y_dot*self.forces.summed_Fy)

        self.forces.Fm = np.sqrt(np.power(self.forces.summed_Fx,2)
                                 +np.power(self.forces.summed_Fy,2)
                                 )
        
        self.forces.tau = self.forces.xpin*self.forces.summed_Fy-self.forces.ypin*self.forces.summed_Fx
        # Calculate the mean normal force on the crank pin
        # This assumes a quasi-steady bearing where the film is well-behaved
        self.forces.mean_Fm = np.trapz(self.forces.Fm[_slice], self.t[_slice])/(2*pi)
        self.forces.mean_Fr = np.trapz(self.forces.Fr[_slice], self.t[_slice])/(2*pi)
        self.forces.mean_Ft = np.trapz(self.forces.Ft[_slice], self.t[_slice])/(2*pi)
        self.forces.mean_tau = np.trapz(self.forces.tau[_slice], self.t[_slice])/(2*pi)
                
        #: The interial forces on the orbiting scroll
        self.forces.inertial = self.orbiting_scroll_mass * self.omega**2 * self.geo.ro / 1000
        
        self.forces.corr_Fm = self.forces.mean_Fm + self.forces.inertial
        
        print 'self.forces.mean_Ft', self.forces.mean_Ft
        print 'self.forces.mean_Fr', self.forces.mean_Fr
        print 'self.forces.mean_Fz', self.forces.mean_Fz
        print 'self.forces.corr_Fz', self.forces.corr_Fz
        
        
        
    def IsentropicNozzleFM(self,*args,**kwargs):
        """
        A thin wrapper around the base class function for pickling purposes
        """
        return PDSimCore.IsentropicNozzleFM(self,*args,**kwargs)
        
        