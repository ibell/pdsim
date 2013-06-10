from __future__ import division
##--- package imports
from PDSim.misc.datatypes import arraym
from PDSim.core.containers import ControlVolume
from PDSim.flow.flow import FlowPath
from PDSim.core.core import PDSimCore
from PDSim.flow import flow_models
from PDSim.plot.plots import debug_plots
from PDSim.core.bearings import journal_bearing
from PDSim.scroll import scroll_geo
from _scroll import _Scroll

##--- non-package imports
import warnings
from scipy.optimize import fsolve, newton
from CoolProp.CoolProp import Props
from CoolProp import State
from math import pi,cos
import numpy as np
import copy
import types

class struct(object):
    pass

class Scroll(PDSimCore, _Scroll):
    """
    This is a python class that implements functionality for a scroll compressor
    
    It is inherited from the PDSimCore class
    """
    
    # Bind some of the cython methods from the base class so that they can
    # pickle properly
    # To add method to existing instance, see http://stackoverflow.com/questions/972/adding-a-method-to-an-existing-object
    #RadialLeakage = _Scroll.RadialLeakage
    
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
        return scroll_geo.SA(theta,self.geo,Vremove = self.geo.Vremove)[0:2]
        
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
            A :class:`State <CoolProp.State.State>` instance for the inlet to the scroll set.  Can be approximate
        outletState
            A :class:`State <CoolProp.State.State>` instance for the outlet to the scroll set.  Can be approximate
            
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
    
    def auto_add_leakage(self,flankFunc,radialFunc, radialFunc_kwargs = {}):
        """
        Add all the leakage terms for the compressor
        
        Parameters
        ----------
        flankFunc : function
            The function to be used for the flank leakage path
        radialFunc : function
            The function to be used for the radial leakage path
        radialFunc_kwargs : dict
            Dictionary of terms to be passed to the radial leakage function
            
        """
        
        #Do the flank leakages
        self.auto_add_flank_leakage(flankFunc)
        #Do the radial leakages
        self.auto_add_radial_leakage(radialFunc, radialFunc_kwargs)
        
    def auto_add_radial_leakage(self, radialFunc, radialFunc_kwargs):
        """
        A function to add all the radial leakage terms
        
        Parameters
        ----------
        radialFunc : function
            The function that will be called for each radial leakage
        radialFunc_kwargs : dict
            Dictionary of terms to be passed to the radial leakage function
        """
        #Get all the radial leakage pairs
        pairs = scroll_geo.radial_leakage_pairs(self.geo)
        
        #Loop over all the radial leakage pairs possible for the given geometry
        for pair in pairs:
            self.add_flow(FlowPath(key1=pair[0],
                                   key2=pair[1],
                                   MdotFcn=radialFunc,
                                   MdotFcn_kwargs = radialFunc_kwargs
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
    
    def calculate_scroll_mass(self):
        """
        Calculate the mass and centroid of the orbiting scroll
        
        Returns
        -------
        mtotal : float
            Total mass of the orbiting scroll (including any additional mass added at centroid)
            
        zcm__thrust_surface : float
            The location of the centroid above the thrust surface. z = 0 is at the thrust surface, and positive z direction is towards the orbiting scroll 
        """
        tplate = getattr(self.mech,'scroll_plate_thickness')
        rho = getattr(self.mech,'scroll_density')
        mplus = getattr(self.mech,'scroll_added_mass')
        Lbearing = getattr(self.mech,'L_crank_bearing')
        Dijournal = getattr(self.mech,'D_crank_bearing')
        Dplate = getattr(self.mech,'scroll_plate_diameter')
        Dojournal = 1.5*Dijournal
        
        Vwrap,cx,cy = scroll_geo.scroll_wrap(self.geo)
        
        mwrap = rho * Vwrap
        mplate = rho * pi * tplate * Dplate**2/4.0
        mjournal = rho * pi * Lbearing * (Dojournal**2-Dijournal**2)/4.0
        mtotal = mwrap + mplate + mjournal + mplus
        
        zwrap = tplate + self.geo.h/2
        zplate = tplate/2
        zjournal = -Lbearing/2
        
        zcm__thrust_surface = (mwrap*zwrap + mjournal*zjournal + mplate*zplate)/mtotal
        
        return mtotal, zcm__thrust_surface

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
        This function evaluates the anti-derivative of the differential of wall heat transfer, and returns the amount of scroll-wall heat transfer in kW
        
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
    
    def heat_transfer_callback(self, theta):
        """
        The scroll simulation heat transfer callback for HT to the fluid in the 
        chambers
        
        ``heat_transfer_callback`` for ``PDSimCore.derivs`` must be of the 
        form::
        
            heat_transfer_callback(theta)
            
        but we need to get the inlet and outlet states to get the linear 
        temperature profile in the scroll wrap. Thus we wrap the callback 
        we would like to call in this function that allows us to determine
        the inlet and outlet state at run-time.
        """
        State_inlet = self.Tubes.Nodes[self.key_inlet]
        State_outlet = self.Tubes.Nodes[self.key_outlet]
        return self._heat_transfer_callback(theta, State_inlet, State_outlet)
    
    def _heat_transfer_callback(self, theta, State_inlet, State_outlet, HTC_tune = 1.0):
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
        
        Q = []
        for key in self.CVs.exists_keys:
            Q.append(self.calcHT(theta,key,HTC_tune,dT_dphi,phim))
        return arraym(Q)
        
    def step_callback(self,t,h,Itheta):
        """
        Here we test whether the control volumes need to be
        a) Merged
        b) Adjusted because you are at the discharge angle
        
        """ 
        #This gets called at every step, or partial step
        self.theta=t
        
        def angle_difference(angle1,angle2):
            # Due to the periodicity of angles, you need to handle the case where the
            # angles wrap around - suppose theta_d is 6.28 and you are at an angles of 0.1 rad
            #, the difference should be around 0.1, not -6.27
            # 
            # This brilliant method is from http://blog.lexique-du-net.com/index.php?post/Calculate-the-real-difference-between-two-angles-keeping-the-sign
            # and the comment of user tk
            return (angle1-angle2+pi)%(2*pi)-pi
        
        def IsAtMerge(eps = 0.001, eps_d1_higher=0.01,eps_dd_higher=0.00001):
            pressures = [self.CVs['d1'].State.p,
                         self.CVs['d2'].State.p,
                         self.CVs['dd'].State.p]
            p_max = max(pressures)
            p_min = min(pressures)
            if abs(p_min/p_max-1)<eps_dd_higher:
                return True
            # For over compression cases, the derivatives don't tend to drive
            # the pressures together, and as a result you need to relax the 
            # convergence quite a bit
            elif angle_difference(t, scroll_geo.theta_d(self.geo))>1.2 and abs(p_min/p_max-1)<eps_d1_higher:
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
            #print 'At the discharge angle'
            ########################
            #Reassign chambers
            ########################
            #Find chambers with a discharge_becomes flag
            for key in self.CVs.exists_keys:
                if self.CVs[key].discharge_becomes in self.CVs.keys:
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
            # That is why it is 2.2e-8 rather than 2.0e-8
            h=2.2e-10
            disable=True
       
        elif self.CVs['d1'].exists and IsAtMerge():
            
            #Build the volume vector using the old set of control volumes (pre-merge)
            V,dV=self.CVs.volumes(t)
            
            print 'merging'
            
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
                
                self.update_existence()
                
                #Re-calculate the CV
                V,dV=self.CVs.volumes(t)
                self.T[self.CVs.exists_indices,Itheta] = self.CVs.T
                self.p[self.CVs.exists_indices,Itheta] = self.CVs.p
                self.m[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)*V
                self.rho[self.CVs.exists_indices,Itheta] = arraym(self.CVs.rho)
                
            else:
                raise NotImplementedError('no flooding yet')
            disable=True 
              
        elif t>self.theta_d:
            self.__before_discharge1__=False
            disable=False
            
        return disable,h
        
    def crank_bearing(self, W):
        
        JB = journal_bearing(r_b = self.mech.D_crank_bearing/2,
                             L = self.mech.L_crank_bearing,
                             omega = self.omega,
                             W = W,
                             c = self.mech.c_crank_bearing,
                             eta_0 = self.mech.mu_oil
                             )
        self.losses.crank_bearing_dict = JB
    
        return JB['Wdot_loss']/1000.0
        
    def upper_bearing(self, W):
        """
        Moment balance around the upper bearing gives the force for
        the lower bearing.  Torques need to balance around the upper bearing
        """
        
        JB = journal_bearing(r_b = self.mech.D_upper_bearing/2,
                             L = self.mech.L_upper_bearing,
                             omega = self.omega,
                             W = W,
                             c = self.mech.c_upper_bearing,
                             eta_0 = self.mech.mu_oil
                             )
        self.losses.upper_bearing_dict = JB

        return JB['Wdot_loss']/1000.0
    
    def lower_bearing(self, W):
        """
        Moment balance around the upper bearing gives the force for
        the lower bearing.  Torques need to balance around the upper bearing
        """
        
        JB = journal_bearing(r_b = self.mech.D_lower_bearing/2,
                             L = self.mech.L_lower_bearing,
                             omega = self.omega,
                             W = W,
                             c = self.mech.c_lower_bearing,
                             eta_0 = self.mech.mu_oil
                             )
        self.losses.lower_bearing_dict = JB

        return JB['Wdot_loss']/1000.0
    
    def thrust_bearing(self):
        """
        The thrust bearing analysis
        """
        from PDSim.core.bearings import thrust_bearing
        V = self.geo.ro*self.omega
        #Use the corrected force to account for the decrease in back area due to the bearing
        N = self.forces.summed_Fz*1000 #[N]
        TB = thrust_bearing(mu = self.mech.thrust_friction_coefficient,
                            V = V,
                            N = N)
        self.losses.thrust_bearing_dict = TB
        return TB['Wdot_loss']/1000.0
    
    def mechanical_losses(self, shell_pressure = 'low:shell'):
        """
        Calculate the mechanical losses in the bearings
        
        Parameters
        ----------
            shell_pressure : string, 'low', 'low:shell', 'mid', or 'high'
            
            low uses the upstream pressure of the machine,
            
            low:shell uses the pressure after the inlet tube
            
            mid uses the average of upstream and downstream pressures
            
            high uses the pressure downstream of the machine 

        """
        
        #inlet pressure [kPa]
        inlet_pressure = self.Tubes.Nodes[self.key_inlet].p
        outlet_pressure = self.Tubes.Nodes[self.key_outlet].p
        
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
        else:
            raise KeyError("keyword argument shell_pressure must be one of 'low', 'low:shell', 'mid' or 'high'; received '"+str(shell_pressure)+"'")
        
        #Calculate the force terms: force profiles, mean values, etc. 
        self.calculate_force_terms(orbiting_back_pressure = back_pressure)
        
        if not hasattr(self,'losses'):
            self.losses = struct()
            
        if not hasattr(self.mech,'journal_tune_factor'):
            warnings.warn('mech.journal_tune_factor not provided; using 1.0')
            self.mech.journal_tune_factor = 1.0
            
        if hasattr(self.mech,'detailed_analysis') and self.mech.detailed_analysis == True:
            self.detailed_mechanical_analysis()
            return self.forces.Wdot_total_mean
        else:
            #Conduct the calculations for the bearings [N]
            W_OSB = np.sqrt((self.forces.summed_Fr + self.forces.inertial)**2+self.forces.summed_Ft**2)*1000
            
            _slice = range(self.Itheta+1)
            
            self.losses.crank_bearing = np.zeros_like(W_OSB)
            self.losses.upper_bearing = np.zeros_like(W_OSB)
            self.losses.lower_bearing = np.zeros_like(W_OSB)
            
            for i in _slice:
                self.losses.crank_bearing[i] = self.crank_bearing(W_OSB[i])
                self.losses.upper_bearing[i] = self.upper_bearing(W_OSB[i]*(1+1/self.mech.L_ratio_bearings))
                self.losses.lower_bearing[i] = self.lower_bearing(W_OSB[i]/self.mech.L_ratio_bearings)
                
            self.losses.thrust_bearing = self.thrust_bearing()
            
            theta = self.t[_slice]
            theta_range = theta[-1]-theta[0]
            
            # Sum up each loss v. theta curve
            self.losses.summed = self.losses.crank_bearing + self.losses.upper_bearing + self.losses.lower_bearing + self.losses.thrust_bearing
            
            # Get the mean losses over one cycle
            self.losses.bearings  = np.trapz(self.losses.summed[_slice], theta)/theta_range
            
            print 'mechanical losses: ', self.losses.bearings
            return self.losses.bearings #[kW]
    
    def post_cycle(self):
        #Run the base-class method to set HT terms, etc.
        PDSimCore.post_cycle(self)
        #Calculate the mechanical and motor losses 
        self.lump_energy_balance_callback()
        #Update the heat transfer to the gas in the shell
        self.suction_heating()
        
    def ambient_heat_transfer(self, Tshell):
        """
        The amount of heat transfer from the compressor to the ambient
        """
        return self.h_shell*self.A_shell*(Tshell-self.Tamb)
    
    def initial_motor_losses(self, eta_a = 0.8):
        """
        Assume an adiabatic efficiency to estimate the motor power and motor losses
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
        
        # Get an initial guess before running at all for the motor losses.
        self.initial_motor_losses() #set the parameter self.motor.losses
        
        # Run the suction heating code
        self.suction_heating()
        
        # Calculate the dischare port free area
        #self.cache_discharge_port_blockage()
        
        # Call the base class function        
        PDSimCore.pre_run(self)
        
    def cache_discharge_port_blockage(self, xport = None, yport = None):
        """
        Precalculate the discharge port blockage using the polygon math module
        
        This is computationally quite expensive, which is why it is done only once
        and then interpolated within.
        """
        
        import matplotlib.pyplot as plt
        from PDSim.misc.polymath import Polygon, PolygonOperator
        
        if xport is None and yport is None:
            xport = self.geo.xa_arc1 + 0.9*self.geo.ra_arc1*np.cos(np.linspace(0,2*pi,100))
            yport = self.geo.ya_arc1 + 0.9*self.geo.ra_arc1*np.sin(np.linspace(0,2*pi,100))
        
        pport = Polygon(xport,yport)
        
        xinvi, yinvi = scroll_geo.coords_inv(np.linspace(self.geo.phi_is+pi/2,self.geo.phi_is,25), self.geo, 0, flag="fi")
        xarc1 = self.geo.xa_arc1 + self.geo.ra_arc1*np.cos(np.linspace(self.geo.t2_arc1,self.geo.t1_arc1,75))
        yarc1 = self.geo.ya_arc1 + self.geo.ra_arc1*np.sin(np.linspace(self.geo.t2_arc1,self.geo.t1_arc1,75))
        xinvo, yinvo = scroll_geo.coords_inv(np.linspace(self.geo.phi_os,self.geo.phi_os+3*pi/2,25), self.geo, 0, flag="fo")
        
        # TODO: Add line and arc2 if they exist
        t,A=[],[]
        for theta in np.linspace(0,2*pi,100):
            
            THETA = self.geo.phi_ie-pi/2.0-theta
        
            xscroll = -np.r_[xinvi,xarc1,xinvo]+self.geo.ro*np.cos(THETA)
            yscroll = -np.r_[yinvi,yarc1,yinvo]+self.geo.ro*np.sin(THETA)
            
            pscroll = Polygon(xscroll, yscroll)
            
#            plt.plot(xport, yport)
#            plt.plot(xscroll,yscroll)
#            
            s = 0
            for poly in PolygonOperator(pport, pscroll).Only1():
                s += poly.area()
                
            t.append(theta)
            A.append(s)
            
        plt.plot(t,A)
        plt.show()
        
    def guess_lump_temps(self, T0):
        """
        Guess the temperature of the lump
        
        Parameters
        ----------
        T0 : float
            First guess for temperature [K]
        """
        
        # First try to just alter the lump temperature with the gas heat transfer
        # rate fixed
        
        def OBJECTIVE(x):
            self.Tlumps[0] = x
            #Run the tubes
            for tube in self.Tubes:
                tube.TubeFcn(tube)
            #
            return self.lump_energy_balance_callback()[0]
        
        print OBJECTIVE(T0-50)
        print OBJECTIVE(T0+50)
        return newton(OBJECTIVE,T0)
        
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
        self.mech.Wdot_losses = self.mechanical_losses('low:shell') 
        
        # Heat transfer with the ambient; Qamb is positive if heat is being removed, thus flip the sign
        Qnet -= self.Qamb
        
        Qnet += self.mech.Wdot_losses
        # Heat transfer with the gas in the working chambers.  mean_Q is positive
        # if heat is transfered to the gas in the working chamber, so flip the 
        # sign for the lump
        Qnet -= self.HTProcessed.mean_Q
        
        #Shaft power from forces on the orbiting scroll from the gas in the pockets [kW]
        self.Wdot_forces = self.omega*self.forces.mean_tau
        
        self.Wdot_mechanical = self.Wdot_pv + self.mech.Wdot_losses
        
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
            raise AttributeError('motor.type must be one of "const_eta_motor" or "motor_map"')

#        print 'mean_Q', self.HTProcessed.mean_Q
#        print 'self.forces.mean_Fm', self.forces.mean_Fm
#        print 'self.forces.inertial', self.forces.inertial
#        print 'self.Qamb', self.Qamb
#        print 'self.Wdot_forces', self.Wdot_forces
#        print 'self.Wdot_pv', self.Wdot_pv 
#        print 'self.losses.bearings', self.losses.bearings
#        print 'self.Wdot_mechanical', self.Wdot_mechanical
        
        #Motor losses [kW]
        self.motor.losses = self.Wdot_mechanical*(1/self.eta_motor-1)
        
        #Electrical Power
        self.Wdot_electrical = self.Wdot_mechanical + self.motor.losses
        
        if hasattr(self,'Wdot_i'):
            #Overall isentropic efficiency
            self.eta_oi = self.Wdot_i/self.Wdot_electrical
        
#        #Set the heat input to the suction line
#        self.suction_heating()
        
        print 'At this iteration'
        print '    Electrical power:', self.Wdot_electrical,'kW'
        print '    Mass flow rate:', self.mdot,'kg/s'
        if hasattr(self,'Wdot_i'):
            print '    Over. isentropic:', self.eta_oi,'-'
        if hasattr(self,'eta_v'):
            print '    Volumetric:', self.eta_v,'-'
        
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
     
    def SA_S1(self, FlowPath, X_d=1.0,**kwargs):
        """
        A wrapper for the flow between the suction area and the S1 chamber
        
        Notes
        -----
        If geo.phi_ie_offset is greater than 0, the offset geometry will be 
        used to calculate the flow area.  Otherwise the conventional analysis 
        will be used.
        """
        if abs(self.geo.phi_ie_offset) > 1e-12:
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
        
    def Injection_to_Comp(self,FlowPath,phi,inner_outer,check_valve = False, A = 7e-6, **kwargs):
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
        
        """
        #1. Figure out what CV is connected to the port
        partner_key = self._get_injection_CVkey(phi, self.theta, inner_outer)
        FlowPath.A = A
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
        # then there is no back flow, and hence no flow at all
        elif check_valve:
            if FlowPath.key_down.startswith('inj'):
                return 0.0
            
#                #This will be negative
#                DELTAp = FlowPath.State_down.p - FlowPath.State_up.p
#            else:
#                #This will be positive
#                DELTAp = FlowPath.State_up.p - FlowPath.State_down.p 
#            
#            # Using an approximation to a Heaviside step function to close off the
#            # port gradually and improve numerical convergence due to continuous
#            # first derivative
#            if -10 < DELTAp < 10.0:
#                FlowPath.A *=  1/(1+np.exp(-10*(DELTAp-2)))
#            elif DELTAp < -10.0:
#                return 0.0
        
        mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                            FlowPath.State_up,
                                            FlowPath.State_down)
        return mdot
        
    def calculate_force_terms(self,
                              orbiting_back_pressure=None):
        """
        Calculate the force profiles, mean forces, moments, etc.
        
        Parameters
        ----------
        orbiting_back_pressure : float, or class instance
            If a class instance, must provide a function __call__ that takes as its first input the Scroll class
        
        """
        if not isinstance(orbiting_back_pressure,float):
            raise ValueError('calculate_force_terms must get a float back pressure for now')
            
        if not hasattr(self.mech,'scroll_plate_thickness'):
            warnings.warn('"mech.scroll_plate_thickness" not found, using 2*scroll wrap thickness')
            self.mech.scroll_plate_thickness = 2*self.geo.t
        
        if not hasattr(self.mech,'scroll_zcm__thrust_surface'):
            warnings.warn('"mech.scroll_zcm__thrust_surface" not found, using 0')
            self.mech.scroll_zcm__thrust_surface = 0
        
        self.forces = struct()
        
        #Get the slice of indices that are in use.  At the end of the simulation
        #execution this will be the full range of the indices, but when used
        # at intermediate iterations it will be a subset of the indices
        _slice = range(self.Itheta+1)
        
        t = self.t[_slice]
        
        ####################################################
        #############  Inertial force components ###########
        ####################################################
        
        #: The magnitude of the inertial forces on the orbiting scroll [kN]
        self.forces.inertial = self.mech.orbiting_scroll_mass * self.omega**2 * self.geo.ro / 1000
        
        ####################################################
        #############  Normal force components #############
        ####################################################
        # The force of the gas in each chamber pushes the orbiting scroll away
        # from the working chambers
        # It is only the active slice
        self.forces.Fz = self.p[:,_slice]*self.V[:,_slice]/self.geo.h
        
        #Remove all the NAN placeholders and replace them with zero values
        self.forces.Fz[np.isnan(self.forces.Fz)] = 0
        #Sum the terms for the applied gas force from each of the control volumes
        self.forces.summed_Fz = np.sum(self.forces.Fz, axis = 0) #kN
#        
        # The back gas pressure on the orbiting scroll pushes the scroll back down
        # Subtract the back pressure from all the control volumes 
        self.forces.Fbackpressure = orbiting_back_pressure*self.V[:,_slice]/self.geo.h
        
        #Remove all the NAN placeholders and replace them with zero values
        self.forces.Fbackpressure[np.isnan(self.forces.Fbackpressure)] = 0
        self.forces.summed_Fbackpressure = np.sum(self.forces.Fbackpressure, axis = 0) #kN
        
        self.forces.summed_Fz -= self.forces.summed_Fbackpressure
        
        # Add the net axial force generated by the gas at the top of the scroll wrap
        self.forces.Faxial_involute = self.scroll_involute_axial_force(t) 
        
        self.forces.summed_Faxial_involute = np.sum(self.forces.Faxial_involute, axis = 0)
        
        self.forces.summed_Fz += self.forces.summed_Faxial_involute
        
        # Calculate the mean axial force
        self.forces.mean_Fz = np.trapz(self.forces.summed_Fz, t)/(2*pi)
        
        ####################################################
        #############  "Radial" force components ###########
        ####################################################
        self.forces.Fx = np.zeros((self.CVs.N,len(self.t[_slice])))
        self.forces.Fy = np.zeros_like(self.forces.Fx)
        self.forces.fxp = np.zeros_like(self.forces.Fx)
        self.forces.fyp = np.zeros_like(self.forces.Fx)
        self.forces.cx = np.zeros_like(self.forces.Fx)
        self.forces.cy = np.zeros_like(self.forces.Fx)
        self.forces.Mz = np.zeros_like(self.forces.Fx)
        
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
        for CVkey in self.CVs.keys:
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
                self.forces.fxp[I,:] = [comp['fx_p'] for comp in geo_components]
                self.forces.fyp[I,:] = [comp['fy_p'] for comp in geo_components]
                self.forces.Fx[I,:] = [comp['fx_p'] for comp in geo_components]*p
                self.forces.Fy[I,:] = [comp['fy_p'] for comp in geo_components]*p
                self.forces.cx[I,:] = [comp['cx'] for comp in geo_components]
                self.forces.cy[I,:] = [comp['cy'] for comp in geo_components]
                self.forces.Mz[I,:] = [comp['M_O_p'] for comp in geo_components]*p
        
        # The thrust load from JUST the working chambers
        self.forces.summed_gas_Fz = np.sum(self.forces.Fz, axis = 0)
        
        # Point of action of the thrust forces - weighted by the axial force
        self.forces.cx_thrust = np.sum(self.forces.cx*self.forces.Fz, axis = 0) / self.forces.summed_gas_Fz
        self.forces.cy_thrust = np.sum(self.forces.cy*self.forces.Fz, axis = 0) / self.forces.summed_gas_Fz
        
        self.forces.THETA = self.geo.phi_ie-pi/2-self.t[_slice]
        # Position of the pin as a function of crank angle
        self.forces.xpin = self.geo.ro*np.cos(self.forces.THETA)
        self.forces.ypin = self.geo.ro*np.sin(self.forces.THETA)
        
        # Remove all the NAN placeholders
        self.forces.Fx[np.isnan(self.forces.Fx)] = 0
        self.forces.Fy[np.isnan(self.forces.Fy)] = 0
        self.forces.Fz[np.isnan(self.forces.Fz)] = 0
        self.forces.Mz[np.isnan(self.forces.Mz)] = 0
        
        # Sum the terms at each crank angle
        self.forces.summed_Fx = np.sum(self.forces.Fx, axis = 0) #kN
        self.forces.summed_Fy = np.sum(self.forces.Fy, axis = 0) #kN
        
        # Moment around the x axis and y-axis from the thrust load caused by the working chambers
        # Sign convention on force is backwards here (in this case, forces pointing down are positive
        # r is vector from (xpin,ypin) to (cx,cy), F is -Fz*k
        #        |   i          j      k   |  
        #   M =  | cx-xpin   cy-ypin   0   | = - (cy-ypin) * Fz *i + (cx- xpin) * Fz * j
        #        |   0          0     -Fz  |
        #   
        self.forces.Mx = -(self.forces.cy - self.forces.ypin)*self.forces.Fz
        self.forces.My = +(self.forces.cx - self.forces.xpin)*self.forces.Fz
        
        # Moment around the x-axis and y-axis from the applied gas load on the orbiting scroll wrap relative to the thrust plane
        #        |   i    j     k  |  
        #   M =  |   0    0     z  | = - Fy*z *i + Fx * z * j
        #        |   Fx   Fy    0  |
        #
        self.forces.Mx += -self.forces.Fy*(self.mech.scroll_plate_thickness + self.geo.h/2)
        self.forces.My += self.forces.Fx*(self.mech.scroll_plate_thickness + self.geo.h/2)
        
        self.forces.summed_Mz = np.sum(self.forces.Mz, axis = 0) #kN-m
        self.forces.summed_Mx = np.sum(self.forces.Mx, axis = 0) #kN-m
        self.forces.summed_My = np.sum(self.forces.My, axis = 0) #kN-m
        
        # Magnitude of the overturning moment generated by the gas forces (Fr, Ftheta, Fz)
        self.forces.Moverturn_gas = np.sqrt(self.forces.summed_Mx**2+self.forces.summed_My**2)
        
        # The moments from the backpressure acting on the back side of the orbiting scroll
        # They must be added on separately because otherwise they are added in NCV times,
        # which is not the proper behavior
        # r is vector from (xpin,ypin) to (0,0), F is Fbp*k
        #        |   i        j      k   |  
        #   M =  | 0-xpin   0-ypin   0   | = -ypin * Fbp *i + xpin * Fbp * j
        #        |   0        0     Fbp  |
        #        
        self.forces.summed_Mx += -self.forces.ypin*self.forces.summed_Fbackpressure
        self.forces.summed_My +=  self.forces.xpin*self.forces.summed_Fbackpressure
#        
        # Moment around the x-axis and y-axis from the inertial force of the orbiting scroll
        # They must be added on separately because otherwise they are added in NCV times,
        # which is not the proper behavior
        #
        # Components of the inertial force in the x- and y-directions [kN]
        Fcx = self.forces.inertial*np.cos(self.forces.THETA)
        Fcy = self.forces.inertial*np.sin(self.forces.THETA)
        
        # Inertial overturning moment acts through the center of mass of the orbiting scroll
        # Moment around the x-axis and y-axis from the applied gas load on the orbiting scroll wrap relative to the thrust plane
        #        |   i      j     k  |  
        #   M =  |   0      0     z  | = - Fcy*z *i + Fcx * z * j
        #        |   Fcx   Fcy    0  |
        #
        self.forces.summed_Mx += -Fcy*(self.mech.scroll_zcm__thrust_surface)
        self.forces.summed_My += Fcx*(self.mech.scroll_zcm__thrust_surface)
        
        # Moment around the x-axis and y-axis from the journal bearing applied forces on the orbiting scroll
        self.forces.Fx_bearing_simple = -(np.cos(self.forces.THETA)*self.forces.inertial + self.forces.summed_Fx)
        self.forces.Fy_bearing_simple = -(np.sin(self.forces.THETA)*self.forces.inertial + self.forces.summed_Fy)
        #        |   i      j     k   |  
        #   M =  |   0      0     -h  | = Fy*h *i - Fx * h * j
        #        |   Fx     Fy    0   |
        #
        self.forces.summed_Mx +=  self.forces.Fy_bearing_simple*(self.mech.D_crank_bearing/2)
        self.forces.summed_My += -self.forces.Fx_bearing_simple*(self.mech.D_crank_bearing/2)
        
        self.forces.Moverturn_thrust = np.sum(self.forces.Fz*np.sqrt((self.forces.cy - self.forces.ypin)**2+(self.forces.cx - self.forces.xpin)**2),axis=0)
        self.forces.Moverturn_gas = (self.mech.scroll_plate_thickness + self.geo.h/2)*np.sum(np.sqrt(self.forces.Fx**2+self.forces.Fy**2),axis=0)
        self.forces.Moverturn_inertia = (self.mech.scroll_zcm__thrust_surface)*np.sqrt(Fcy**2+Fcx**2)
        self.forces.Moverturn_bearing = (self.mech.D_crank_bearing/2)*np.sqrt(self.forces.Fx_bearing_simple**2+self.forces.Fy_bearing_simple**2)
        self.forces.Moverturn = np.sqrt(self.forces.summed_Mx**2+self.forces.summed_My**2)
        
#        import matplotlib.pyplot as plt
#        plt.plot(t,self.forces.Moverturn_thrust,'b')
#        plt.plot(t,self.forces.Moverturn_gas,'r')
#        plt.plot(t,self.forces.Moverturn_inertia,'g')
#        plt.plot(t,self.forces.Moverturn_bearing,'k')
#        plt.show()
        
        # Center of reaction in the global coordinate system
        self.forces.x_thrust_reaction = -self.forces.summed_My/self.forces.summed_Fz+self.forces.xpin #Fz is positive if pushing UP on the OS
        self.forces.y_thrust_reaction = -self.forces.summed_Mx/self.forces.summed_Fz+self.forces.ypin #Fz is positive if pushing UP on the OS
        self.forces.r_thrust_reaction = np.sqrt(self.forces.x_thrust_reaction**2+ self.forces.y_thrust_reaction**2)
        
        #Calculate the radial force on the crank pin at each crank angle
        #The radial component magnitude is just the projection of the force onto a vector going from origin to center of orbiting scroll
        self.forces.Fr = (np.cos(self.forces.THETA)*self.forces.Fx + np.sin(self.forces.THETA)*self.forces.Fy)
        #
        #Components of the unit vector in the direction of rotation
        x_dot = +np.sin(self.forces.THETA)
        y_dot = -np.cos(self.forces.THETA)
        # Direction of rotation is opposite the positive theta direction, so need to flip sign for Ft
        self.forces.Ft = -(x_dot*self.forces.Fx+y_dot*self.forces.Fy)
        
        #Remove all the NAN placeholders
        self.forces.Fr[np.isnan(self.forces.Fr)] = 0
        self.forces.Ft[np.isnan(self.forces.Ft)] = 0
        #Sum the terms at each crank angle
        self.forces.summed_Fr = np.sum(self.forces.Fr,axis = 0) #kN
        self.forces.summed_Ft = np.sum(self.forces.Ft,axis = 0) #kN

        self.forces.Fm = np.sqrt(self.forces.summed_Fx**2+self.forces.summed_Fy**2)
        
        self.forces.tau = self.forces.xpin*self.forces.summed_Fy-self.forces.ypin*self.forces.summed_Fx
        # Calculate the mean normal force on the crank pin
        # This assumes a quasi-steady bearing where the film is well-behaved
        self.forces.mean_Fm = np.trapz(self.forces.Fm, self.t[_slice])/(2*pi)
        self.forces.mean_Fr = np.trapz(self.forces.summed_Fr, self.t[_slice])/(2*pi)
        self.forces.mean_Ft = np.trapz(self.forces.summed_Ft, self.t[_slice])/(2*pi)
        self.forces.mean_tau = np.trapz(self.forces.tau, self.t[_slice])/(2*pi)
        self.forces.mean_Mz = np.trapz(self.forces.summed_Mz, self.t[_slice])/(2*pi)
        
    def detailed_mechanical_analysis(self):
        """
        In this function the detailed mechanical analsysis is carried out.
        
        Notes
        -----
        This function implements the method of the documentation on the orbiting scroll
        forces
        
        The applied force on each of the bearings can be given from::
        
                        **|  |**    -----
                        **|  |**     |
                      |      |       | z_upper
                      |      |       |
                    **|      |**     |
                    **|      |**    -----   
                    **|      |**     |
                      |      |       |
                      |      |       |
                      |      |       |
                      |      |       | z_lower
                      |      |       |
                      |      |       |
                      |      |       |
                    **|      |**     |
                    **|      |**   -----
        
        and the ratio of bearing distances is given by ``mech.L_ratio_bearings`` which is z_lower/z_upper
        """
        
        if not hasattr(self,'losses'):
            self.losses = struct()
            
        muthrust = self.mech.thrust_friction_coefficient
        beta = self.mech.oldham_rotation_beta
        mu1 = mu2 = mu3 = mu4 = mu5 = self.mech.oldham_key_friction_coefficient
        r1 = r2 = r3 = r4 = self.mech.oldham_ring_radius
        w1 = w2 = w3 = w4 = self.mech.oldham_key_width
        mOR = self.mech.oldham_mass
        mOS = self.mech.orbiting_scroll_mass
        wOR = self.mech.oldham_thickness
        hkeyOR = self.mech.oldham_key_height
        
        # Fill in terms if they are not provided for backwards compatability
        for term in ['pin1_ybeta_offset','pin2_ybeta_offset','pin3_xbeta_offset','pin4_xbeta_offset']:
            if not hasattr(self.mech, term):
                warnings.warn('"mech.'+term+'" not found, using 0')
                setattr(self.mech,term,0)

        F1_ybeta_offset = self.mech.pin1_ybeta_offset
        F2_ybeta_offset = self.mech.pin2_ybeta_offset
        F3_xbeta_offset = self.mech.pin3_xbeta_offset
        F4_xbeta_offset = self.mech.pin4_xbeta_offset
        
        # Gravitional acceleration
        g = 9.80665 
        
        _slice = range(self.Itheta+1)
        theta = self.t[_slice]
        
        # The initial guesses for the moment generated by the journal bearing - 
        # it should be positive since Ms is negative and Ms and M_B act in 
        # opposite directions
        self.forces.F_B0 = np.sqrt((self.forces.summed_Fr + self.forces.inertial)**2+self.forces.summed_Ft**2)
        self.forces.mu_B = np.zeros_like(self.forces.F_B0)
        self.forces.mu_Bupper = np.zeros_like(self.forces.F_B0)
        self.forces.mu_Blower = np.zeros_like(self.forces.F_B0)
        self.forces.M_B0 = self.forces.mu_B*self.mech.D_crank_bearing/2*self.forces.F_B0
        
        THETA = self.geo.phi_ie-pi/2-theta
        vOR_xbeta = self.geo.ro*self.omega*(np.sin(THETA)*np.cos(beta)-np.cos(THETA)*np.sin(beta)) #Velocity of the oldham ring in the xbeta direction
        aOR_xbeta = self.geo.ro*self.omega**2*(-np.cos(THETA)*np.cos(beta)-np.sin(THETA)*np.sin(beta))
        UPSILON = vOR_xbeta/np.abs(vOR_xbeta)
        vOS_ybeta = -self.geo.ro*self.omega*(np.sin(THETA)*np.sin(beta)+np.cos(THETA)*np.cos(beta)) #Velocity of the orbiting scroll in the y_beta direction
        PSI = vOS_ybeta/np.abs(vOS_ybeta)
        vOS = self.geo.ro*self.omega
        aOS_x = -self.geo.ro*self.omega**2*np.cos(THETA)
        aOS_y = -self.geo.ro*self.omega**2*np.sin(THETA)
        
        Nsteps = self.Itheta+1
        A = np.zeros((4,4,Nsteps))
        b = np.zeros((4,Nsteps))
        self.forces.Fkey = np.zeros((4,Nsteps))
        
        # Make a matrix stack where each entry in the third index corresponds to a 4x4 matrix of the terms
        # Oldham x-beta direction
        A[0,0,:] = -mu1*UPSILON
        A[0,1,:] = -mu2*UPSILON
        A[0,2,:] = 1
        A[0,3,:] = -1
        b[0,:] = mOR*aOR_xbeta/1000
        
        # Oldham ybeta direction
        A[1,0,:] = 1    
        A[1,1,:] = -1 
        A[1,2,:] = -mu3*PSI
        A[1,3,:] = -mu4*PSI
        b[1,:] = 0
            
        # Oldham moments around the central z-direction axis
        A[2,0,:] = r1-mu1*UPSILON*(w1/2-F1_ybeta_offset)
        A[2,1,:] = r2+mu2*UPSILON*(w2/2+F2_ybeta_offset)
        A[2,2,:] = -r3+mu3*PSI*(w3/2-F3_xbeta_offset)
        A[2,3,:] = -r4-mu4*PSI*(w4/2+F4_xbeta_offset)
        b[2,:] = 0
        
        # Orbiting scroll moments around the central axis
        A[3,0,:] = 0
        A[3,0,:] = 0
        A[3,0,:] = r3-mu3*PSI*(w3/2-F3_xbeta_offset)
        A[3,0,:] = r4+mu4*PSI*(w4/2+F4_xbeta_offset)
        
        # Use the initial guess for the bearing moments and applied force
        self.forces.M_B = self.forces.M_B0
        self.forces.F_B = self.forces.F_B0
        
        step = 1
        
        # In the first step, we use the guess values for the bearing normal force, 
        # in the second step we use the back-calculated bearing normal force 
        while step <= 2:
            
            #Calculate the friction coefficient for each bearing
            for i in _slice:
                self.crank_bearing(self.forces.F_B[i]*1000)
                self.forces.mu_B[i] = self.losses.crank_bearing_dict['f']
                self.upper_bearing(self.forces.F_B[i]*1000*(1+1/self.mech.L_ratio_bearings))
                self.forces.mu_Bupper[i] = self.losses.upper_bearing_dict['f']
                self.lower_bearing(self.forces.F_B[i]*1000/self.mech.L_ratio_bearings)
                self.forces.mu_Blower[i] = self.losses.lower_bearing_dict['f']
                
            self.forces.M_B = self.forces.mu_B*self.mech.D_crank_bearing/2*self.forces.F_B
            
            # This term depends on M_B which is re-calculated at each iteration.  All other terms are independent of M_B
            b[3,:] = -self.forces.summed_Mz-self.forces.M_B-muthrust*(self.forces.summed_My*np.cos(THETA)-self.forces.summed_Mx*np.sin(THETA))
            
            # Walk through the stack of arrays, and obtain F1,F2,F3,F4 for each
            # crank angle
            for i in _slice:
                self.forces.Fkey[:,i] = np.linalg.solve(A[:,:,i], b[:,i])
        
                #debug_plots(self)
            F1, F2, F3, F4 = self.forces.Fkey
        
            # Bearing forces on the scroll re-calculated based on force balances in the x- and y-axes
            Fbx = mOS*aOS_x/1000-muthrust*self.forces.summed_Fz*np.sin(THETA)+mu3*PSI*F3*np.sin(beta)+mu4*PSI*F4*np.sin(beta)-F4*np.cos(beta)+F3*np.cos(beta)-self.forces.summed_Fx
            Fby = mOS*aOS_y/1000-muthrust*self.forces.summed_Fz*np.cos(THETA)-mu3*PSI*F3*np.cos(beta)-mu4*PSI*F4*np.cos(beta)-F4*np.sin(beta)+F3*np.sin(beta)-self.forces.summed_Fy
            Fbold = np.sqrt(Fbx**2+Fby**2)
            self.forces.M_B = self.forces.mu_B*self.mech.D_crank_bearing/2*Fbold
            self.forces.M_Bupper = self.forces.mu_Bupper*self.mech.D_upper_bearing/2*Fbold*(1.0+1.0/self.mech.L_ratio_bearings)
            self.forces.M_Blower = self.forces.mu_Blower*self.mech.D_lower_bearing/2*Fbold/self.mech.L_ratio_bearings
            
            self.forces.F_B = Fbold
            
            step += 1
        
        self.forces.Wdot_F1 = np.abs(F1*vOR_xbeta*mu1)
        self.forces.Wdot_F2 = np.abs(F2*vOR_xbeta*mu2)
        self.forces.Wdot_F3 = np.abs(F3*vOS_ybeta*mu3)
        self.forces.Wdot_F4 = np.abs(F4*vOS_ybeta*mu4)
        self.forces.Wdot_OS_journal = np.abs(self.omega*self.forces.M_B)*self.mech.journal_tune_factor
        self.forces.Wdot_upper_journal = np.abs(self.omega*self.forces.M_Bupper)*self.mech.journal_tune_factor
        self.forces.Wdot_lower_journal = np.abs(self.omega*self.forces.M_Blower)*self.mech.journal_tune_factor
        
        self.forces.Wdot_thrust = np.abs(muthrust*self.forces.summed_Fz*vOS)
        
        self.forces.Wdot_total = (  self.forces.Wdot_F1
                                  + self.forces.Wdot_F2
                                  + self.forces.Wdot_F3
                                  + self.forces.Wdot_F4
                                  + self.forces.Wdot_OS_journal
                                  + self.forces.Wdot_upper_journal 
                                  + self.forces.Wdot_lower_journal 
                                  + self.forces.Wdot_thrust)
        
        self.forces.Wdot_total_mean = np.trapz(self.forces.Wdot_total, theta)/(2*pi)
        print self.forces.Wdot_total_mean,'average mechanical losses'
            
        import matplotlib.pyplot as plt
            
#        fig = plt.figure()
#        fig.add_subplot(111)
#        plt.plot(theta,self.forces.Wdot_F1,label='F1')
#        plt.plot(theta,self.forces.Wdot_F2,label='F2')
#        plt.plot(theta,self.forces.Wdot_F3,label='F3')
#        plt.plot(theta,self.forces.Wdot_F4,label='F4')
#        plt.plot(theta,self.forces.Wdot_OS_journal,label='journal')
#        plt.plot(theta,self.forces.Wdot_thrust,label='thrust')
#        plt.plot(theta,self.forces.Wdot_total,lw=3)
#        plt.legend()
#        plt.show()
#        
#        plt.close('all')
#        fig = plt.figure()
#        fig.add_subplot(141)
#        plt.gca().set_title('Oldham acceleration [g]')
#        plt.plot(theta,(F[2,:].T-F[3,:].T)*1000/9.81/mOR,'o',lw=2)
#        plt.plot(theta, aOR_xbeta/9.81,'-',lw=2)
#        fig.add_subplot(142)
#        plt.plot(theta,F.T)
#        plt.gca().set_title('Key forces [kN]')
#        fig.add_subplot(143)
#        plt.plot(theta,self.forces.summed_Ft)
#        plt.gca().set_title('Tangential force [kN]')
#        fig.add_subplot(144)
#        plt.plot(theta,self.forces.summed_Fr)
#        plt.gca().set_title('Radial force [kN]')
#        plt.show()
        
    def scroll_involute_axial_force(self, theta, p_backpressure = 0):
        """
        Calculate the axial force generated by the pressure distribution 
        along the top of the scroll wrap 
        
        Pressure along inner and outer walls is considered to act over one-half 
        of the thickness of the scroll wrap.
        
        Notes
        -----
        The following assumptions are employed:
        
        1. Involute extended to the base circle to account for discharge region
        2. Half of the width of the scroll is assumed to see the upstream pressure and the other half sees the downstream pressure
        
        The length of an involute section can be given by
        
        .. math::
            
            s = r_b\\left(\\frac{\\phi_2^2-\\phi_1^2}{2}-\\phi_{0}(\\phi_2-\\phi_1)\\right)
        
        Returns
        -------
        F: numpy array
            Axial force matrix from the top of the scroll generated by each control volume [kN]
        """
        
        _slice = range(len(theta))
        
        #Get the break angle (simplified solution)
        phi_s_sa = self.geo.phi_oe-pi
        
        #Get the number of compression chambers in existence at each crank angle
        nC = np.array([scroll_geo.getNc(t,self.geo) for t in theta])
        
        F = np.zeros_like(self.p)
        F = F[:,_slice]
        
        #Parameters for the SA chamber
        phi2 = self.geo.phi_oe
        phi1 = phi_s_sa
        ds_SA = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_o0*(phi2-phi1))
        I = self.CVs.index('sa')
        F[I,:] = ds_SA*self.geo.t/2*(self.p[I,_slice]-p_backpressure)
        
        #Parameters for the S1 chamber
        phi2 = phi_s_sa
        phi1 = phi_s_sa-theta
        ds_S1 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_o0*(phi2-phi1))
        I = self.CVs.index('s1')
        F[I,:] = ds_S1*self.geo.t/2*(self.p[I,_slice]-p_backpressure)
        
        # Parameters for the S2 chamber
        phi2 = self.geo.phi_ie
        phi1 = self.geo.phi_ie-theta
        ds_S2 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_i0*(phi2-phi1))
        I = self.CVs.index('s2')
        F[I,:] = ds_S2*self.geo.t/2*(self.p[I,_slice]-p_backpressure)
        
        for I in range(1, scroll_geo.nC_Max(self.geo)+1):
            phi2 = self.geo.phi_oe-pi-theta-2*pi*(I-1)
            phi1 = self.geo.phi_oe-pi-theta-2*pi*(I)
            ds_C1 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_o0*(phi2-phi1))
            ICV = self.CVs.index('c1.'+str(I))
            F[ICV,:] = ds_C1*self.geo.t/2*(self.p[ICV, _slice]-p_backpressure)
            
            phi2 = self.geo.phi_ie-theta-2*pi*(I-1)
            phi1 = self.geo.phi_ie-theta-2*pi*(I)
            ds_C2 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_i0*(phi2-phi1))
            ICV = self.CVs.index('c2.'+str(I))
            F[ICV,:] = ds_C2*self.geo.t/2*(self.p[ICV, _slice]-p_backpressure)
        
        phi2 = self.geo.phi_oe-pi-theta-2*pi*(nC)
        phi1 = self.geo.phi_os
        ds_D1 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_o0*(phi2-phi1))
        ICV = self.CVs.index('d1')
        F[ICV,:] = ds_D1*self.geo.t/2*(self.p[ICV,_slice]-p_backpressure)
        
        phi2 = self.geo.phi_ie-theta-2*pi*(nC)
        phi1 = self.geo.phi_is
        ds_D2 = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_i0*(phi2-phi1))
        ICV = self.CVs.index('d2')
        F[ICV,:] = ds_D2*self.geo.t/2*(self.p[ICV,_slice]-p_backpressure)
        
        phi2 = self.geo.phi_is
        phi1 = self.geo.phi_i0
        ds_DD = self.geo.rb*(0.5*(phi2**2-phi1**2)-self.geo.phi_i0*(phi2-phi1))
        ICV = self.CVs.index('dd')
        F[ICV,:] = ds_DD*self.geo.t/2*(self.p[ICV,_slice]-p_backpressure)
        
        ICV = self.CVs.index('ddd')
        F[ICV,:] = (ds_D1+ds_D2+ds_DD)*self.geo.t/2*(self.p[ICV,_slice]-p_backpressure)
        
        # Remove all the nan placeholders
        F[np.isnan(F)] = 0
        
        return F
        
    def IsentropicNozzleFMSafe(self,*args,**kwargs):
        """
        A thin wrapper around the base class function for pickling purposes
        """
        return PDSimCore.IsentropicNozzleFMSafe(self, *args, **kwargs)
        
    def IsentropicNozzleFM(self,*args,**kwargs):
        """
        A thin wrapper around the base class function for pickling purposes
        """
        return PDSimCore.IsentropicNozzleFM(self, *args, **kwargs)
    
    def attach_HDF5_annotations(self, fName):
        """
        In this function, annotations can be attached to each HDF5 field
        
        Here we add other scroll-specific terms
        
        Parameters
        ----------
        fName : string
            The file name for the HDF5 file that is to be used
        """ 
        #Use the base annotations
        PDSimCore.attach_HDF5_annotations(self, fName)
        
        attrs_dict = {
                '/forces/F_B':'The normal force applied to the orbiting scroll journal bearing [kN]',
                '/forces/F_B0':'The normal force applied to the orbiting scroll journal bearing [kN]',
                '/forces/Fkey':'The forces applied at each key of the Oldham ring [kN]',
                '/forces/Fbackpressure':'The force applied to the orbiting scroll due to the back pressure [kN]',
                '/forces/Fm':'The normal force applied to the orbiting scroll journal bearing [kN]',
                '/forces/Fr':'The radial gas force on the orbiting scroll [kN]',
                '/forces/Ft':'The tangential gas force on the orbiting scroll [kN]',
                '/forces/Fx':'The gas force on the orbiting scroll in the x-direction [kN]',
                '/forces/Fy':'The gas force on the orbiting scroll in the y-direction [kN]',
                '/forces/Fz':'The gas force on the orbiting scroll in the negative z-direction [kN]',
                '/forces/M_B':'The journal bearing moment on the orbiting scroll in the positive x-direction [kN-m]',
                '/forces/M_B0':'The journal bearing moment on the orbiting scroll in the positive x-direction [kN-m]',
                '/forces/Moverturn':'The magnitude of the overturning moment on the orbiting scroll [kN-m]',
                '/forces/Mx':'The overturning moment around the x-axis [kN-m]',
                '/forces/My':'The overturning moment around the y-axis [kN-m]',
                '/forces/Mz':'The spinning moment from the gas around the z-axis [kN-m]',
                '/forces/THETA':'The shifted angle that is used to locate the center of the orbiting scroll [rad]',
                '/forces/Wdot_F1':'The frictional dissipation at key 1 of Oldham ring [kW]',
                '/forces/Wdot_F2':'The frictional dissipation at key 2 of Oldham ring [kW]',
                '/forces/Wdot_F3':'The frictional dissipation at key 3 of Oldham ring [kW]',
                '/forces/Wdot_F4':'The frictional dissipation at key 4 of Oldham ring [kW]',
                '/forces/Wdot_OS_journal':'The frictional dissipation at journal bearing of orbiting scroll [kW]',
                '/forces/Wdot_upper_journal':'The frictional dissipation at upper journal bearing [kW]',
                '/forces/Wdot_lower_journal':'The frictional dissipation at lower journal bearing [kW]',
                '/forces/Wdot_thrust':'The frictional dissipation from thrust bearing [kW]',
                '/forces/Wdot_total':'The frictional dissipation of bearings and friction [kW]',
                '/forces/Wdot_total_mean':'The mean frictional dissipation of bearings and friction over one rotation[kW]',
                '/forces/cx':'The x-coordinate of the centroid of each control volume [m]',
                '/forces/cx_thrust':'Effective x-coordinate of the centroid of all chambers [m]',
                '/forces/cy':'The y-coordinate of the centroid of each control volume [m]',
                '/forces/cy_thrust':'Effective y-coordinate of the centroid of all chambers [m]',
                '/forces/fxp':'Fx/p from geometric analysis for each control volume [kN/kPa]',
                '/forces/fyp':'Fy/p from geometric analysis for each control volume [kN/kPa]',
                '/forces/inertial':'Magnitude of inertial force (m*omega^2*r) [kN]',
                '/forces/mean_Fm':'Mean of Fm over one rotation [kN]',
                '/forces/mean_Fr':'Mean of Fr over one rotation [kN]',
                '/forces/mean_Ft':'Mean of Ft over one rotation [kN]',
                '/forces/mean_Fz':'Mean of Fz over one rotation [kN]',
                '/forces/mean_Mz':'Mean of Mz over one rotation [kN-m]',
                '/forces/mean_tau':'Mean of tau over one rotation [kN-m]',
                '/forces/summed_Fr':'Summation of CV contributions to Fr [kN]',
                '/forces/summed_Ft':'Summation of CV contributions to Ft [kN]',
                '/forces/summed_Fx':'Summation of CV contributions to Fx [kN]',
                '/forces/summed_Fy':'Summation of CV contributions to Fy [kN]',
                '/forces/summed_Fz':'Summation of CV contributions to Fz [kN]',
                '/forces/summed_Moverturn':'Summation of CV contributions to overturning moment [kN-m]',
                '/forces/summed_Mx':'Summation of CV contributions to Mx [kN-m]',
                '/forces/summed_My':'Summation of CV contributions to My [kN-m]',
                '/forces/summed_Mz':'Summation of CV contributions to Mz [kN-m]',
                '/forces/summed_gas_Fz':'Summation of CV contributions to Fz (only from the CV) [kN-m]',
                '/forces/tau':'Torque generated by gas around the central axis of shaft [kN-m]',
                '/forces/xpin':'x-coordinate of the orbiting scroll center [m]',
                '/forces/ypin':'y-coordinate of the orbiting scroll center [m]',
                '/mech/D_crank_bearing':'Diameter of the crank journal bearing [m]',
                '/mech/D_lower_bearing':'Diameter of the lower journal bearing [m]',
                '/mech/D_upper_bearing':'Diameter of the upper journal bearing [m]',
                '/mech/L_crank_bearing':'Length of the crank journal bearing [m]',
                '/mech/L_lower_bearing':'Length of the lowe journal bearing [m]',
                '/mech/L_ratio_bearings':'Ratio of the distances from the upper bearing to the crank bearing [m]',
                '/mech/L_upper_bearing':'Length of the upper journal bearing [m]',
                '/mech/c_crank_bearing':'Clearance (D/2-rb) of the crank journal bearing [m]',
                '/mech/c_lower_bearing':'Clearance (D/2-rb) of the lower journal bearing [m]',
                '/mech/c_upper_bearing':'Clearance (D/2-rb) of the upper journal bearing [m]',
                '/mech/detailed_analysis':'True if detailed mechanical analysis is being used',
                '/mech/journal_tune_factor':'Tuning factor that muliplies losses in each journal bearing [-]',
                '/mech/mu_oil':'Viscosity of the oil [Pa-s]',
                '/mech/oldham_key_friction_coefficient':'Friction coefficient for all keys in Oldham ring [-]',
                '/mech/oldham_key_height':'Height of each key of Oldham ring [m]',
                '/mech/oldham_key_width':'Width of each key of Oldham ring [m]',
                '/mech/oldham_mass':'Mass of Oldham ring [kg]',
                '/mech/oldham_ring_radius':'Radius of Oldham ring [m]',
                '/mech/oldham_rotation_beta':'Angle between Oldham sliding axis and x-axis [radian]',
                '/mech/oldham_thickness':'Thickness of the main ring of Oldham ring [m]',
                '/mech/orbiting_scroll_mass':'Mass of orbiting scroll [kg]',
                '/mech/scroll_density':'Scroll density [kg]',
                '/mech/scroll_plate_thickness':'Scroll plate thickness [m]',
                '/mech/thrust_ID':'Thrust bearing inner diameter [m]',
                '/mech/thrust_OD':'Thrust bearing outer diameter [m]',
                '/mech/thrust_friction_coefficient':'Thrust bearing friction coefficient [m]'
                }
        
        import h5py
        hf = h5py.File(fName,'a')
        
        for k, v in attrs_dict.iteritems():
            dataset = hf.get(k)
            if dataset is not None:
                dataset.attrs['note'] = v
        hf.close()
        
        