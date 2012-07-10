##--- package imports
from PDSim.core.containers import ControlVolume
from PDSim.flow.flow import FlowPath
from ..core.core import PDSimCore
from PDSim.misc._listmath import listm
import scroll_geo
from PDSim.flow import flow_models
from ..plot.plots import debug_plots

##--- non-package imports
from scipy.optimize import fsolve
from CoolProp.CoolProp import UseSinglePhaseLUT
from CoolProp import State
from math import pi,cos
from _scroll import _Scroll
import copy

class Scroll(PDSimCore):
    """
    This is a Python class that extends functionality in the GeometryCore.cpp SWIGGed C++ code
    
    It is inherited from the PDModelCore class
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
        #py_dict.update(_Scroll.__cdict__(self))
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
        return (3*pi-2*self.geo.phi_ie+self.geo.phi_i0+self.geo.phi_o0)/(-2*self.geo.phi_os-3*pi+self.geo.phi_i0+self.geo.phi_o0)
    
    def V_injection(self, theta, V_tube = None):
        """
        Volume of injection tube
        
        The tube volume can either be given by the keyword argument V_tube 
        (so you can easily have more than one injection tube), or it can be 
        provided by setting the attribute V_inj_tube and NOT providing the V_tube
        argument
        """
        if V_tube is None:
            return self.V_inj_tube, 0.0
        else:
            return V_tube, 0.0
        
    def V_sa(self,theta,full_output=False):
        """
        Wrapper around the Cython code for sa calcs
        
        theta: angle in range [0,2*pi]
        
        Returns
        -------
        
        """
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return scroll_geo.SA(theta,self.geo)[0:2],HTangles
        else:
            return scroll_geo.SA(theta,self.geo)[0:2]
        
    def V_s1(self,theta,full_output=False):
        """
        Wrapper around the Cython code for Vs1_calcs
        
        theta: angle in range [0,2*pi]
        """
        if full_output==True:
            HTangles = {'1_i':self.geo.phi_ie,
                        '2_i':self.geo.phi_ie-theta,
                        '1_o':scroll_geo.phi_s_sa(theta,self.geo),
                        '2_o':self.geo.phi_oe-pi-theta}
            return scroll_geo.S1(theta,self.geo)[0:2],HTangles
        else:
            return scroll_geo.S1(theta,self.geo)[0:2]
        
    def V_s2(self,theta,full_output=False):
        """
        Wrapper around the Cython code for Vs1_calcs
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        if full_output==True:
            HTangles = {'1_i':self.geo.phi_ie,
                        '2_i':self.geo.phi_ie-theta,
                        '1_o':scroll_geo.phi_s_sa(theta,self.geo),
                        '2_o':self.geo.phi_oe-pi-theta}
            return scroll_geo.S2(theta,self.geo)[0:2],HTangles
        else:
            return scroll_geo.S2(theta,self.geo)[0:2]
    
    def V_c1(self,theta,alpha=1,full_output=False):
        """
        Wrapper around the Cython code for C1
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return scroll_geo.C1(theta,alpha,self.geo)[0:2],HTangles
        else:
            return scroll_geo.C1(theta,alpha,self.geo)[0:2]
        
    def V_c2(self,theta,alpha=1,full_output=False):
        """
        Wrapper around the Cython code for C2
        
        theta: angle in range [0,2*pi]
        alpha: index of compression chamber pair; 1 is for outermost set
        """
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return scroll_geo.C2(theta,alpha,self.geo)[0:2],HTangles
        else:
            return scroll_geo.C2(theta,alpha,self.geo)[0:2]
        
    def V_d1(self,theta,full_output=False):
        """
        Wrapper around the compiled code for D1
        
        theta: angle in range [0,2*pi]
        """
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return scroll_geo.D1(theta,self.geo)[0:2],HTangles
        else:
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
        if full_output==True:
            HTangles = {'1_i':None,'2_i':None,'1_o':None,'2_o':None}
            return scroll_geo.D2(theta,self.geo)[0:2],HTangles
        else:
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
        
    def set_scroll_geo(self,Vdisp,Vratio,Thickness,OrbitingRadius,phi_i0=0.0,phi_os=0.3):
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
        """
        
        phi_i0=0
        phi_is=pi
        phi_os=0.3
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
                p1 = inletState.p
                T1 = inletState.T
                k = inletState.cp/inletState.cv
                V1 = self.V_s1(2*pi)[0]
                V2 = self.V_c1(0,alpha)[0]
                p2 = p1 * (V1/V2)**k
                T2 = T1 * (V1/V2)**(k-1)
                initState=State.State(inletState.Fluid,dict(T=T2,P=p2)).copy()
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
    
    def auto_add_leakage(self,flankFunc,radialFunc,):
        
        # First do the flank leakage terms since they are easier to handle
        # Always a s1-c1 leakage and s2-c2 leakage
        self.add_flow(FlowPath(key1='s1',key2='c1.1',MdotFcn=flankFunc))
        self.add_flow(FlowPath(key1='s2',key2='c2.1',MdotFcn=flankFunc))
        
        # Only add the DDD-S1 and DDD-S2 flow path if there is one set of
        # compression chambers.  You can provide a different flank leakage
        # function if desired by setting the 
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
                self.add_flow(FlowPath(key1=keyc1,key2='c1.'+str(alpha+1),MdotFcn=flankFunc))
                self.add_flow(FlowPath(key1=keyc2,key2='c2.'+str(alpha+1),MdotFcn=flankFunc))
                
            elif alpha==nCmax:
                #Leakage between the discharge region and the innermost chamber
                self.add_flow(FlowPath(key1=keyc1,key2='ddd',MdotFcn=flankFunc))
                self.add_flow(FlowPath(key1=keyc2,key2='ddd',MdotFcn=flankFunc))
    
#    def RadialLeak(self, geo, phi_max, phi_min):
#        flowVec->A[*count]=(geo->delta_radial)*(geo->rb)*(1.0/2.0*(phi_max*phi_max-phi_min*phi_min)-(geo->phi.phi_fi0)*(phi_max-phi_min))
    
    def heat_transfer_callback(self,theta,**kwargs):
        
#        def func(CV):
#            #Join them all into one tuple
#            V_dV,HTangles=CV.V_dV(theta,full_output=True,**CV.V_dV_kwargs)
#            return V_dV+(HTangles,)
#        #Loop over the control volumes that exist 
#        V_dV_HTangles=map(func,self.CVs.exists_CV)
#        V,dV,HTangles=zip(*V_dV_HTangles)
#        def calcHT(V_dV_HTangles):
#            V=V_dV_HTangles[0]
#            HTangles=V_dV_HTangles[2]
#            if None in HTangles.values():
#                return 0.0
#            phi_1_i=HTangles['1_i']
#            phi_2_i=HTangles['2_i']
#            phi_1_o=HTangles['1_o']
#            phi_2_o=HTangles['2_o']
#            
#            A_plate = V/self.geo.h
#            A_wall_in = self.geo.h * self.geo.rb * ( (phi_1_i**2-phi_2_i**2)/2.0 - self.geo.phi_i0*(phi_1_i-phi_2_i) )
#            A_wall_out = self.geo.h * self.geo.rb * ( (phi_1_o**2 - phi_2_o**2)/2.0 -self.geo.phi_o0*(phi_1_o-phi_2_o) )
#            return 0.0
#        Q=map(calcHT,V_dV_HTangles)
        return listm([0.0]*self.CVs.Nexist)
        
    def step_callback(self,t,h,Itheta,**kwargs):
        """
        Here we test whether the control volumes need to be
        a) Merged
        b) Adjusted because you are at the discharge angle
        
        """ 
        #This gets called at every step, or partial step
        self.theta=t
        
        def IsAtMerge(eps_d1_higher=0.005,eps_dd_higher=0.0001):
            if self.CVs['d1'].State.p>self.CVs['dd'].State.p and abs(self.CVs['d1'].State.p/self.CVs['dd'].State.p-1)<eps_d1_higher:
                print 'Merged with d1 higher'
                return True
            elif self.CVs['d1'].State.p<self.CVs['dd'].State.p and abs(self.CVs['d1'].State.p/self.CVs['dd'].State.p-1)<eps_dd_higher:
                print 'Merged with dd higher'
                return True
            else:
                return False
            
        disable=False
        
        if t<self.theta_d<t+h and self.__before_discharge2__==False:
            #Take a step almost up to the discharge angle
            disable=True
            h=self.theta_d-t-1e-8
            self.__before_discharge2__=True
        elif self.__before_discharge2__==True:
            #At the discharge angle
            print 'At the discharge angle'
            #-----------------
            #Reassign chambers
            #-----------------
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
            h=2.2e-8
            disable=True
       
        elif self.CVs['d1'].exists and IsAtMerge():
            
            #Build the volume vector using the old set of control volumes (pre-merge)
            V,dV=self.CVs.volumes(t)
            
            if self.__hasLiquid__==False:

                #Density
                rhod1=self.CVs['d1'].State.rho
                rhod2=self.CVs['d2'].State.rho
                rhodd=self.CVs['dd'].State.rho
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
                u_before=ud1*rhod1*Vd1+ud2*rhod2*Vd2+udd*rhodd*Vdd
                rhoddd=m/Vddd
                T=(Td1*Vd1+Td2*Vd2+Tdd*Vdd)/Vddd
                self.CVs['ddd'].State.update({'T':T,'D':rhoddd})
                u_after=self.CVs['ddd'].State.u*self.CVs['ddd'].State.rho*Vddd
                
                DeltaU=m*(u_before-u_after)
                if abs(DeltaU)>1e-6:
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
        
    def mechanical_losses(self):
        print 'temporary mechanical losses'
        return 0.2 #[kW]
    
    def ambient_heat_transfer(self,Tshell):
        """
        The amount of heat transfer from the compressor to the ambient
        """
        return 0.001*(Tshell-self.Tamb)
    
    def lump_energy_balance_callback(self):
        
        #For the single lump
        # terms are positive if heat transfer is TO the lump
        Qnet=0.0
        for Tube in self.Tubes:
            Qnet-=Tube.Q
        Qnet+=self.mechanical_losses()-self.ambient_heat_transfer(self.Tlumps[0])
        return [Qnet]
    
    def TubeCode(self,Tube,**kwargs):
        Tube.Q = flow_models.IsothermalWallTube(Tube.mdot,Tube.State1,Tube.State2,Tube.fixed,Tube.L,Tube.ID,T_wall=self.Tlumps[0])
        
    def InjectionTubeFM(self,FlowPath,**kwargs):
        FlowPath.A=pi*0.005**2/4.0
        try:
            return flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
        except ZeroDivisionError:
            return 0.0
    
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
        
    def Inlet_sa(self, FlowPath, X_d=1.0, **kwargs):
        FlowPath.A=X_d*pi*0.03**2/4.0
        return flow_models.IsentropicNozzle(FlowPath.A,
                                            FlowPath.State_up,
                                            FlowPath.State_down)
        
    def Discharge(self, FlowPath, **kwargs):
        FlowPath.A=pi*0.01**2/4.0
        try:
            return flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
        except ZeroDivisionError:
            return 0.0
        
    def SA_S(self, FlowPath, X_d=1.0,**kwargs):
        FlowPath.A=X_d*scroll_geo.Area_s_sa(self.theta, self.geo)
        try:
            mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
        
    def FlankLeakage(self,FlowPath,**kwargs):
        """
        This function 
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
#            return flow_models.IsentropicNozzle(FlowPath.A,
#                                                FlowPath.State_up,
#                                                FlowPath.State_down)
        except ZeroDivisionError:
            return 0.0
        
    def _get_injection_CVkey(self,phi,theta,inner_outer):
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
        
    def Injection_to_Comp(self,FlowPath,phi,inner_outer,**kwargs):
        """
        Function to calculate flow rate with
        
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
        #3. Find the distance of the scroll from the point on the involute
        #   where the port is tangent
        FlowPath.A = pi*0.005**2/4.0
        try:
            mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
        
    def RadialLeakage(self,FlowPath,**kwargs):
        assert (1==0)