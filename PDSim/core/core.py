##-- Non-package imports  --
from __future__ import division
from CoolProp.CoolProp import Props
from CoolProp.State import State
import numpy as np
from math import pi
import textwrap
from scipy.integrate import trapz, simps
from scipy.optimize import newton
import copy
import pylab
from time import clock

##--  Package imports  --
from PDSim.flow._sumterms import setcol, getcol
from PDSim.flow import _flow
from PDSim.flow.flow import FlowPathCollection, FlowPath
from containers import ControlVolumeCollection
from PDSim.plot.plots import debug_plots
from PDSim.misc._listmath import listm 
from PDSim.misc.solvers import Broyden,MultiDimNewtRaph

#An empty class for storage
class struct():
    pass    
        
class TubeCollection(list):
    
    def _Nodes(self):
        """
        Nodes is a dictionary of flow states for any tubes that exist
        """
        list1=[(Tube.key1,Tube.State1) for Tube in self if Tube.exists==True]
        list2=[(Tube.key2,Tube.State2) for Tube in self if Tube.exists==True]
        return dict(list1+list2)
    Nodes=property(_Nodes)
    
class Tube():
    """
    A tube is a component of the model that allows for heat transfer and pressure drop.
    
    With this class, the state of at least one of the points is fixed.  For instance, at the inlet of the compressor, the state well upstream is quasi-steady.
    """
    def __init__(self,key1,key2,L,ID,OD=None,State1=None,State2=None,fixed=-1,TubeFcn=None,mdot=-1,exists=True):
        self.key1=key1
        self.key2=key2
        self.fixed=fixed
        self.exists=exists
        if fixed<0:
            raise AttributeError(textwrap.dedent("""You must provide an integer 
            value for fixed, either 1 for Node 1 fixed, or 2 for Node 2 fixed.  
            You provided None (or didn\'t include the parameter"""))
        if fixed==1 and isinstance(State1,State) and State2==None:
            #Everything good
            self.State1=State1
            self.State2=State(self.State1.Fluid,{'T':self.State1.T,'D':self.State1.rho})
        elif fixed==2 and isinstance(State2,State) and State1==None:
            #Everything good
            self.State2=State2
            self.State1=State(self.State2.Fluid,{'T':self.State2.T,'D':self.State2.rho})
        else:
            raise AttributeError('Incompatibility between the value for fixed and the states provided')
            
        self.TubeFcn=TubeFcn
        if mdot<0:
            self.mdot=0.010
            print('Warning: mdot not provided to Tube class contructor, guess value of '+str(self.mdot)+' kg/s used')
        else:
            self.mdot=mdot
        self.L=L
        self.ID=ID
        self.OD=OD

    
class PDSimCore(object):
    """
    This is the main driver class for the model
    
    This class is not intended to be run on its own.  It must be subclassed and extended to provide functions for mass flow, etc. 
    """
    def __init__(self,stateVariables=None):
        """
        Initialization of the PD Core
        """
        #Initialize the containers to be empty
        self.Valves=[]
        self.CVs=ControlVolumeCollection()
        self.Flows=FlowPathCollection()
        self.FlowStorage=[]
        self.Tubes=TubeCollection()
        self.Tlumps=np.zeros((1,1))
        self.steps=[]
        self.__hasValves__=False
        if isinstance(stateVariables,list):
            self.stateVariables=stateVariables
        else:
            self.stateVariables=['T','D']
        self._want_abort = False
    
    def __get_from_matrices(self,i):
        """
        Get values back from the matrices
        """
        
        if self.__hasLiquid__==True:
            raise NotImplementedError
            #return np.hstack([self.T[:,i],self.m[:,i],self.xL[:,i]])
        else:
            VarList=[]
            exists_indices = self.CVs.exists_indices
            for s in self.stateVariables:
                if s=='T':
                    #VarList+=list(self.T[exists_indices,i])
                    VarList+=getcol(self.T,i,exists_indices)
                elif s=='D':
                    #VarList+=list(self.rho[exists_indices,i])
                    VarList+=getcol(self.rho,i,exists_indices)
                elif s=='M':
                    #VarList+=list(self.rho[exists_indices,i]*self.V_)
                    VarList+=list(getcol(self.rho,i,exists_indices)*self.V_)
                else:
                    raise KeyError
            if self.__hasValves__:
                VarList+=list(self.xValves[:,i])
            return listm(VarList)
        
    def __put_to_matrices(self,x,i):
        """
        Take a stacked set of T,rho,xL or T,rho and put back in numpy arrays
        """
        exists_indices=self.CVs.exists_indices
        Nexist = self.CVs.Nexist
        Ns = len(self.stateVariables)
        if self.__hasLiquid__==True:
            raise ValueError
            self.T[:,i]=x[0:self.NCV]
            self.m[:,i]=x[self.NCV:2*self.NCV]
            self.xL[:,i]=x[2*self.NCV:3*self.NCV]
        else: #(self.__hasLiquid__==False)
            for iS,s in enumerate(self.stateVariables):
                x_=listm(x[iS*self.CVs.Nexist:self.CVs.Nexist*(iS+1)])
                if s=='T':
                    setcol(self.T, i, exists_indices, x_)
                elif s=='D':
                    setcol(self.rho, i, exists_indices, x_)
                elif s=='M':
                    setcol(self.m, i, exists_indices, x_)
            #Left over terms are for the valves
            if self.__hasValves__:
                setcol(self.xValves, i, range(len(self.Valves)*2), listm(x[Ns*Nexist::]))
        
        setcol(self.V, i, exists_indices, self.V_)
        setcol(self.dV, i, exists_indices, self.dV_)
        setcol(self.p, i, exists_indices, self.p_)
        if hasattr(self,'h_'):
            setcol(self.h, i, exists_indices, self.h_)
        setcol(self.Q, i, exists_indices, self.Q_)
    
    def __postprocess_flows(self):
        """
        In this private method, the flows from each of the flow nodes are summed for 
        each step of the revolution, and then averaged flow rates are calculated.
        """
        
        def sum_flows(key,Flows):
            """
            Sum all the terms for a given flow key.  
            
            Flows "into" the node are positive, flows out of the 
            node are negative
            
            Use the code in the Cython module
            """
            return _flow.sum_flows(key, Flows)
            
        def collect_keys(Tubes,Flows):
            """
            Get all the keys for a given collection of flow elements
            """
            keys=[]
            for Tube in Tubes:
                if Tube.key1 not in keys:
                    keys.append(Tube.key1)
                if Tube.key2 not in keys:
                    keys.append(Tube.key2)
            for Flow in Flows:
                if Flow.key1 not in keys:
                    keys.append(Flow.key1)
                if Flow.key2 not in keys:
                    keys.append(Flow.key2)
            return keys
        
        #Get all the nodes that can exist for tubes and CVs
        keys=collect_keys(self.Tubes,self.Flows)
        
        #Get the instantaneous net flow through each node
        #   and the averaged mass flow rate through each node
        self.FlowsProcessed=struct()
        self.FlowsProcessed.summed_mdot={}
        self.FlowsProcessed.summed_mdoth={}
        self.FlowsProcessed.mean_mdot={}
        self.FlowsProcessed.integrated_mdoth={}
        self.FlowsProcessed.integrated_mdot={}
        self.FlowsProcessed.t=self.t[0:self.Ntheta]

        for key in keys:
            # Empty container numpy arrays
            self.FlowsProcessed.summed_mdot[key]=np.zeros((self.Itheta,))
            self.FlowsProcessed.summed_mdoth[key]=np.zeros((self.Itheta,))
            for i in range(self.Ntheta):
                mdot,mdoth=sum_flows(key,self.FlowStorage[i])
                self.FlowsProcessed.summed_mdot[key][i]=mdot
                self.FlowsProcessed.summed_mdoth[key][i]=mdoth
            
            # All the calculations here should be done in the time domain,
            # rather than crank angle.  So convert angle to time by dividing 
            # by omega, the rotational speed in rad/s.
            trange = self.t[self.Ntheta-1]-self.t[0]
            # integrated_mdoth has units of kJ/rev * f [Hz] --> kJ/s or kW
            self.FlowsProcessed.integrated_mdoth[key]=trapz(self.FlowsProcessed.summed_mdoth[key], 
                                                            self.t[0:self.Ntheta]/self.omega)*self.omega/(2*pi)
            # integrated_mdot has units of kg/rev * f [Hz] --> kg/s
            self.FlowsProcessed.integrated_mdot[key]=trapz(self.FlowsProcessed.summed_mdot[key], 
                                                           self.t[0:self.Ntheta]/self.omega)*self.omega/(2*pi)
            self.FlowsProcessed.mean_mdot[key]=np.mean(self.FlowsProcessed.integrated_mdot[key])
            
            
        # Special-case the tubes.  Only one of the nodes can have flow.  
        #   The other one is invariant because it is quasi-steady.
        for Tube in self.Tubes:
            mdot1 = self.FlowsProcessed.mean_mdot[Tube.key1]
            mdot2 = self.FlowsProcessed.mean_mdot[Tube.key2]
            mdot_i1 = self.FlowsProcessed.integrated_mdot[Tube.key1]
            mdot_i2 = self.FlowsProcessed.integrated_mdot[Tube.key2]
            mdoth_i1 = self.FlowsProcessed.integrated_mdoth[Tube.key1]
            mdoth_i2 = self.FlowsProcessed.integrated_mdoth[Tube.key2]
            #Swap the sign so the sum of the mass flow rates is zero
            self.FlowsProcessed.mean_mdot[Tube.key1]-=mdot2
            self.FlowsProcessed.mean_mdot[Tube.key2]-=mdot1
            self.FlowsProcessed.integrated_mdot[Tube.key1]-=mdot_i2
            self.FlowsProcessed.integrated_mdot[Tube.key2]-=mdot_i1
            self.FlowsProcessed.integrated_mdoth[Tube.key1]-=mdoth_i2
            self.FlowsProcessed.integrated_mdoth[Tube.key2]-=mdoth_i1
            
    def add_flow(self,FlowPath):
        #Add FlowPath instance to the list of flow paths
        self.Flows.append(FlowPath)
        
    def add_CV(self,CV):
        """
        Add a control volume to the model
        
        Parameters
        ----------
        CV : Control Volume instance
            An initialized control volume.  See :class:`PDSim.core.containers.ControlVolume`
            
        """

        if CV.key in self.CVs:
            raise KeyError('Sorry but the key for your Control Volume ['+CV.key+'] is already in use')
        
        #Add the CV to the collection
        self.CVs[CV.key]=CV
        
    def add_tube(self,Tube):
        """
        Add a tube to the model.  Alternatively call PDSimCore.Tubes.append(Tube)
        
        Parameters
        ----------
        Tube : Tube instance
            An initialized Tube.  See :class:`PDSim.core.Core.Tube`
        """
        #Add it to the list
        self.Tubes.append(Tube)
        
    def add_valve(self,Valve):
        """
        Add a valve to the model.  Alternatively call PDSimCore.Valve.append(Tube)
        
        Parameters
        ----------
        Valve : ValveModel instance
            An initialized Tube.  See :class:`PDSim.flow.FlowModels.ValveModel`
        """
        #Add it to the list
        self.Valves.append(Valve)
        self.__hasValves__=True
        
    def __pre_run(self):
        #Build the full numpy arrays for temperature, volume, etc.
        self.t=np.zeros((50000,))
        self.T=np.zeros((self.CVs.N,50000))
        self.T.fill(np.nan)
        self.p=self.T.copy()
        self.h = self.T.copy()
        self.m=self.T.copy()
        self.V=self.T.copy()
        self.dV=self.T.copy()
        self.rho=self.T.copy()
        self.Q=self.T.copy()
        self.xL=self.T.copy()
        self.xValves=np.zeros((2*len(self.Valves),50000))
        
        self.CVs.rebuild_exists()
        
        #Initialize the control volumes
        self.__hasLiquid__=False
        
    def __pre_cycle(self, x0 = None):
        self.temp_vectors_list=[]
        def add_thing(name,item):
            #Actually create the array
            setattr(self,name,item)
            #Add the name of the vector to the list (for easy removal)
            self.temp_vectors_list.append(name)
            
        #Build temporary arrays to avoid constantly creating numpy arrays
        add_thing('T_',listm([0.0]*self.CVs.Nexist))
        add_thing('p_',listm([0.0]*self.CVs.Nexist))
        add_thing('V_',listm([0.0]*self.CVs.Nexist))
        add_thing('dV_',listm([0.0]*self.CVs.Nexist))
        add_thing('rho_',listm([0.0]*self.CVs.Nexist))
        add_thing('Q_',listm([0.0]*self.CVs.Nexist))
        
        self.t.fill(np.nan)
        self.T.fill(np.nan)
        self.p.fill(np.nan)
        self.m.fill(np.nan)
        self.V.fill(np.nan)
        self.dV.fill(np.nan)
        self.rho.fill(np.nan)
        self.Q.fill(np.nan)
        self.xL.fill(np.nan)
        
        self.FlowStorage=[]
        
        #Get the volumes at theta=0
        #Note: needs to occur in this function because V needed to calculate mass a few lines below
        V,dV=self.CVs.volumes(0)
        self.t[0]=0
        
        # if x0 is provided, use its values to initialize the chamber states
        if x0 is None:
            # self.CVs.exists_indices is a list of indices of the CV with the same order of entries
            # as the entries in self.CVs.T
            self.T[self.CVs.exists_indices, 0] = self.CVs.T
            self.p[self.CVs.exists_indices, 0] = self.CVs.p
            self.rho[self.CVs.exists_indices, 0] = self.CVs.rho
        else:
            x0_ = copy.copy(x0)
            #x0 is provided
            if self.__hasValves__:
                x0_.extend([0]*len(self.Valves)*2)
            self.__put_to_matrices(x0_, 0)
        
        # Assume all the valves to be fully closed and stationary at the beginning of cycle
        self.xValves[:,0]=0
        
        self.Tubes_hdict={}
        for Tube in self.Tubes:
            self.Tubes_hdict[Tube.key1]=Tube.State1.get_h()
            self.Tubes_hdict[Tube.key2]=Tube.State2.get_h()
        
    def cycle_SimpleEuler(self,N,x_state,tmin=0,tmax=2*pi,step_callback=None,heat_transfer_callback=None,valves_callback=None):
        """
        The simple Euler PDSim ODE integrator
        
        Parameters
        ----------
        N : integer
            Number of steps
        x_state : listm
            The initial values of the variables (only the state variables)
        tmin : float, optional
            Starting value of the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``]
        tmax : float, optional
            Ending value for the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``] 
        step_callback : function, optional 
            A pointer to a function that is called at the beginning of the step.  This function must be of the form:: 
            
                step_callback(t,h,Itheta)
                
            where ``h`` is the step size that the solver wants to use, ``t`` is the current value of the independent variable, and ``Itheta`` is the index in the container variables.  The return values are ignored, so the same callback can be used as for the ``cycle_RK45`` solver 
                
        heat_transfer_callback : function, optional
            If provided, the heat_transfer_callback function should have a format like::
            
                Q_listm=heat_transfer_callback(t)
            
            It should return a ``listm`` instance with the same length as the number of CV in existence.  The entry in the ``listm`` is positive if the heat transfer is TO the fluid in the CV in order to maintain the sign convention that energy (or mass) input is positive.  Will raise an error otherwise
        
        """
        #Do some initialization - create arrays, etc.
        self.__pre_cycle(x_state)
        
        #Step variables
        t0=tmin
        h=(tmax-tmin)/(N-1)
        
        #Get the beginning of the cycle configured
        x_state.extend([0.0]*len(self.Valves)*2)
        xold = listm(x_state[:])
        self.__put_to_matrices(xold, 0)
        
        for Itheta in range(N-1):
            
            #Call the step callback if provided
            if step_callback is not None:
                CVchanged = False
                disable, h = step_callback(t0, h, Itheta, CVchanged = CVchanged)
                if CVchanged:
                    print 'CV have changed'
                    xold = self.__get_from_matrices(Itheta)
            
            # Step 1: derivatives evaluated at old values of t = t0
            f1 = self.derivs(t0, xold, heat_transfer_callback, valves_callback)
            xnew = xold+h*f1
            
            #Store at the current index (first run at index 0)
            self.t[Itheta] = t0
            self.__put_to_matrices(xold, Itheta)
            self.FlowStorage.append(self.Flows.get_deepcopy())
            
            # Everything from this step is finished, now update for the next
            # step coming
            t0+=h
            xold = xnew
            
        #Run this at the end at index N-1
        self.t[N-1]=2*pi
        self.derivs(tmax, xnew, heat_transfer_callback, valves_callback)
        self.__put_to_matrices(xnew, N-1)
        self.FlowStorage.append(self.Flows.get_deepcopy())
        
        # last index is Itheta, number of steps is Itheta+1
        print 'Itheta steps taken', N
        self.Itheta = N
        self.Ntheta = N
        self.__post_cycle()
        
    def cycle_Heun(self,N,x_state, tmin=0,tmax=2*pi,step_callback=None,heat_transfer_callback=None,valves_callback=None):
        """
        Use the Heun method (modified Euler method)
        
        Parameters
        ----------
        N : integer
            Number of steps
        x_state : listm
            The initial values of the variables (only the state variables)
        tmin : float
            Starting value of the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``]
        tmax : float
            Ending value for the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``] 
        step_callback : function, optional 
            A pointer to a function that is called at the beginning of the step.  This function must be of the form:: 
            
                step_callback(t,h,Itheta)
                
            where ``h`` is the step size that the solver wants to use, ``t`` is the current value of the independent variable, and ``Itheta`` is the index in the container variables.  The return values are ignored, so the same callback can be used as for the ``cycle_RK45`` solver 
                
        heat_transfer_callback : function, optional
            If provided, the heat_transfer_callback function should have a format like::
            
                Q_listm=heat_transfer_callback(t)
            
            It should return a ``listm`` instance with the same length as the number of CV in existence.  The entry in the ``listm`` is positive if the heat transfer is TO the fluid in the CV in order to maintain the sign convention that energy (or mass) input is positive.  Will raise an error otherwise
        
        """
        #Do some initialization - create arrays, etc.
        self.__pre_cycle()
        
        #Start at an index of 0
        Itheta=0
        t0=tmin
        h=(tmax-tmin)/(N-1)
        
         #One call to build the flows at start
        x_state.extend([0.0]*len(self.Valves)*2)
        xold = listm(x_state[:])
        self.theta=t0
        self.derivs(t0, xold, heat_transfer_callback, valves_callback)
        self.FlowStorage.append(self.Flows.get_deepcopy())
        self.__put_to_matrices(xold, 0)
    
        for Itheta in range(N):
            if step_callback!=None:
                step_callback(t0,h,Itheta)
                
            xold=self.__get_from_matrices(Itheta)
                        
            # Step 1: derivatives evaluated at old values
            f1=self.derivs(t0,xold,heat_transfer_callback,valves_callback)
            xtemp=xold+h*f1
            
            #Step 2: Evaluated at predicted step
            f2=self.derivs(t0+h,xtemp,heat_transfer_callback,valves_callback)
            xnew = xold + h/2.0*(f1 + f2)
            
            self.t[Itheta+1]=t0+h
            self.__put_to_matrices(xnew,Itheta+1)
            t0+=h
            Itheta+=1
            self.FlowStorage.append(self.Flows.get_deepcopy())
            
        #Run this at the end
        V,dV=self.CVs.volumes(t0)
        self.CVs.updateStates('T',xnew[0:self.CVs.Nexist],'D',xnew[self.CVs.Nexist:2*self.CVs.Nexist])
        
        # last index is Itheta, number of steps is Itheta+1
        print 'Itheta steps taken',Itheta+1
        self.Itheta = Itheta
        self.Ntheta = N
        self.__post_cycle()
        return
        
    def cycle_RK45(self,x_state,hmin=1e-4,tmin=0,tmax=2.0*pi,eps_allowed=1e-10,step_relax=0.9,
              step_callback=None,heat_transfer_callback=None,valves_callback = None,
              UseCashKarp=True,**kwargs):
        """
        
        This function implements an adaptive Runge-Kutta-Feldberg 4th/5th order
        solver for the system of equations
        
        Parameters
        ----------
        x_state : listm
            The initial values of the variables (only the state variables)
        hmin : float
            Minimum step size, something like 1e-5 usually is good.  Don't make this too big or you may not be able to get a stable solution
        tmin : float
            Starting value of the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``]
        tmax : float
            Ending value for the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``]
        eps_allowed : float
            Maximum absolute error of any CV per step allowed.  Don't make this parameter too big or you may not be able to get a stable solution.  Also don't make it too small because then you are going to run into truncation error.
        step_relax : float, optional
            The relaxation factor that is used in the step resizing algorithm.  Should be less than 1.0; you can play with this parameter to improve the adaptive resizing, but should not be necessary.
        step_callback : function, optional 
            A pointer to a function that is called at the beginning of the step.  This function must be of the form:: 
            
                disableAdaptive,h=step_callback(t,h,Itheta)
                
            where ``h`` is the step size that the adaptive solver wants to use, ``t`` is the current value of the independent variable, and ``Itheta`` is the index in the container variables.  The return value ``disableAdaptive`` is a boolean value that describes whether the adaptive method should be turned off for this step ( ``False`` : use the adaptive method), and ``h`` is the step size you want to use.  If you don't want to disable the adaptive method and use the given step size, just::
                
                return False,h
            
            in your code.
        heat_transfer_callback : function, optional
            If provided, the heat_transfer_callback function should have a format like::
            
                Q_listm=heat_transfer_callback(t)
            
            It should return a ``listm`` instance with the same length as the number of CV in existence.  The entry in the ``listm`` is positive if the heat transfer is TO the fluid in the CV in order to maintain the sign convention that energy (or mass) input is positive.  Will raise an error otherwise
        
        Notes
        -----
        
        Mathematically the adaptive solver can be expressed as::
        
            k1=h*dy(xn                                                                   ,t)
            k2=h*dy(xn+1.0/4.0*k1                                                        ,t+1.0/4.0*h)
            k3=h*dy(xn+3.0/32.0*k1+9.0/32.0*k2                                           ,t+3.0/8.0*h)
            k4=h*dy(xn+1932.0/2197.0*k1-7200.0/2197.0*k2+7296.0/2197.0*k3                ,t+12.0/13.0*h)
            k5=h*dy(xn+439.0/216.0*k1-8.0*k2+3680.0/513.0*k3-845.0/4104.0*k4             ,t+h)
            k6=h*dy(xn-8.0/27.0*k1+2.0*k2-3544.0/2565.0*k3+1859.0/4104.0*k4-11.0/40.0*k5 ,t+1.0/2.0*h)

        where the function dy(y,t) returns a vector of the ODE expressions.
        The new value is calculated from::
        
            xnplus=xn+gamma1*k1+gamma2*k2+gamma3*k3+gamma4*k4+gamma5*k5+gamma6*k6

        In the adaptive solver, the errors for a given step can be calculated from::

            error=1.0/360.0*k1-128.0/4275.0*k3-2197.0/75240.0*k4+1.0/50.0*k5+2.0/55.0*k6

        If the maximum absolute error is above allowed error, the step size is decreased and the step is 
        tried again until the error is below tolerance.  If the error is better than required, the step
        size is increased to minimize the number of steps required.
        
        Before the step is run, a callback the ``step_callback`` method of this class is called.  In the ``step_callback`` callback function you can do anything you want, but you must return 
        """
        
        #Do some initialization - create arrays, etc.
        self.__pre_cycle()
        
        #Start at an index of 0
        Itheta=0
        t0=tmin
        h=hmin
        
        gamma1=16.0/135.0
        gamma2=0.0
        gamma3=6656.0/12825.0
        gamma4=28561.0/56430.0
        gamma5=-9.0/50.0
        gamma6=2.0/55.0
        
         #One call to build the flows at start
        x_state.extend([0.0]*len(self.Valves)*2)
        xold = listm(x_state[:])
        self.theta=t0
        self.derivs(t0, xold, heat_transfer_callback, valves_callback)
        self.FlowStorage.append(self.Flows.get_deepcopy())
        self.__put_to_matrices(xold, 0)
        
        #t is the independent variable here, where t takes on values in the bounded range [tmin,tmax]
        while (t0<tmax):
            
            stepAccepted=False
            while not stepAccepted:
                
                #Reset the flag
                disableAdaptive=False
                
                if t0 + h > tmax:
                    disableAdaptive=True
                    h=tmax - t0
            
                if step_callback is not None and not disableAdaptive:
                    #The user has provided a disabling function for the adaptive method
                    #Call it to figure out whether to use the adaptive method or not
                    #Pass it a copy of the compressor class and the current step size
                    #The function can modify anything in the class to change flags, existence, merge, etc.
                    disableAdaptive,h=step_callback(t0,h,Itheta)
                else:
                    disableAdaptive=False
                    
                if h < hmin and not disableAdaptive:
                    #Step is too small, just use the minimum step size
                    h = 1.0*hmin
                    disableAdaptive = True
                    
                xold=self.__get_from_matrices(Itheta)
                
                    
                if not UseCashKarp:
                    ## Using RKF constants ##
                    
                    # Step 1: derivatives evaluated at old values
                    f1=self.derivs(t0,xold,heat_transfer_callback,valves_callback)
                    xnew1=xold+h*(+1.0/4.0*f1)
                    
                    f2=self.derivs(t0+1.0/4.0*h,xnew1,heat_transfer_callback,valves_callback)
                    xnew2=xold+h*(+3.0/32.0*f1+9.0/32.0*f2)
    
                    f3=self.derivs(t0+3.0/8.0*h,xnew2,heat_transfer_callback,valves_callback)
                    xnew3=xold+h*(+1932.0/2197.0*f1-7200.0/2197.0*f2+7296.0/2197.0*f3)
    
                    f4=self.derivs(t0+12.0/13.0*h,xnew3,heat_transfer_callback,valves_callback)
                    xnew4=xold+h*(+439.0/216.0*f1-8.0*f2+3680.0/513.0*f3-845.0/4104.0*f4)
                    
                    f5=self.derivs(t0+h,xnew4,heat_transfer_callback,valves_callback);
                    xnew5=xold+h*(-8.0/27.0*f1+2.0*f2-3544.0/2565.0*f3+1859.0/4104.0*f4-11.0/40.0*f5)
                    
                    #Updated values at the next step
                    f6=self.derivs(t0+h/2.0,xnew5,heat_transfer_callback,valves_callback)
                    
                    xnew=xold+h*(gamma1*f1 + gamma2*f2 + gamma3*f3 + gamma4*f4 + gamma5*f5 + gamma6*f6)
                        
                    error=h*(1.0/360.0*f1-128.0/4275.0*f3-2197.0/75240.0*f4+1.0/50.0*f5+2.0/55.0*f6)
                else:
                    # Step 1: derivatives evaluated at old values
                    f1=self.derivs(t0,xold,heat_transfer_callback,valves_callback)
                    xnew1=xold+h*(+1.0/5.0*f1)
                    
                    f2=self.derivs(t0+1.0/5.0*h,xnew1,heat_transfer_callback,valves_callback)
                    xnew2=xold+h*(+3.0/40.0*f1+9.0/40.0*f2)
    
                    f3=self.derivs(t0+3.0/10.0*h,xnew2,heat_transfer_callback,valves_callback)
                    xnew3=xold+h*(3.0/10.0*f1-9.0/10.0*f2+6.0/5.0*f3)
    
                    f4=self.derivs(t0+3.0/5.0*h,xnew3,heat_transfer_callback,valves_callback)
                    xnew4=xold+h*(-11.0/54.0*f1+5.0/2.0*f2-70/27.0*f3+35.0/27.0*f4)
                    
                    f5=self.derivs(t0+h,xnew4,heat_transfer_callback,valves_callback);
                    xnew5=xold+h*(1631.0/55296*f1+175.0/512.0*f2+575.0/13824.0*f3+44275.0/110592.0*f4+253.0/4096.0*f5)
                    
                    f6=self.derivs(t0+7.0/8.0*h,xnew5,heat_transfer_callback,valves_callback)
                    
                    #Updated values at the next step using 5-th order
                    xnew=xold+h*(37/378*f1 + 250/621*f3 + 125/594*f4 + 512/1771*f6)
                    
                    error = h*(-277/64512*f1+6925/370944*f3-6925/202752*f4-277.0/14336*f5+277/7084*f6)
                
                max_error=np.sqrt(np.sum(np.power(error,2)))
                
                # If the error is too large, make the step size smaller and try
                # the step again
                if (max_error > eps_allowed):
                    if not disableAdaptive:
                        # Take a smaller step next time, try again on this step
                        # But only if adaptive mode is on
#                        print 'downsizing', h,
                        h *= step_relax*(eps_allowed/max_error)**(0.25)
#                        print h, eps_allowed, max_error
                        stepAccepted=False
                    else:
                        # Accept the step regardless of whether the error 
                        # is too large or not
                        stepAccepted = True
                else:
                    stepAccepted = True
            
            # Step has been accepted (or adaptive disabled), write values back to matrices
            self.t[Itheta+1]=t0+h
            self.CVs.updateStates('T',xnew[0:self.CVs.Nexist],'D',xnew[self.CVs.Nexist:2*self.CVs.Nexist])
            self.p_ = listm(self.CVs.p)
            self.__put_to_matrices(xnew, Itheta+1)
            t0+=h
            Itheta+=1
            self.FlowStorage.append(self.Flows.get_deepcopy())
            
            #Some debugging information
            #print t0,h,hmin
            
            #The error is already below the threshold
            if (max_error < eps_allowed and disableAdaptive == False):
#                print 'upsizing',h,
                #Take a bigger step next time, since eps_allowed>max_error
                h *= step_relax*(eps_allowed/max_error)**(0.2)
#                print h, eps_allowed, max_error
        #Run this at the end
        V,dV=self.CVs.volumes(t0)
        self.CVs.updateStates('T',xnew[0:self.CVs.Nexist],'D',xnew[self.CVs.Nexist:2*self.CVs.Nexist])
        
        # last index is Itheta, number of steps is Itheta+1
        print 'Itheta steps taken', Itheta+1
        self.Itheta=Itheta
        self.Ntheta = Itheta
        self.__post_cycle()
        
    def __post_cycle(self):
        """
        This stuff all happens at the end of the cycle.  It is a private method not meant to be called externally
        """
        for name in self.temp_vectors_list:
            delattr(self,name)
            
        self.Wdot_pv = 0.0
        for CVindex in range(self.p.shape[0]):
            y = self.p[CVindex,0:self.Ntheta]
            x = self.V[CVindex,0:self.Ntheta]
            self.Wdot_pv += -trapz(y, x)*self.omega/(2*pi)
            
        self.__postprocess_flows()
    
    def check_abort(self):
        #If you received an abort request, set a flag
        if self.pipe_abort.poll() and self.pipe_abort.recv():
            print 'received an abort request'
            self._want_abort = True
            
        return self._want_abort
        
    def solve(self,
              key_inlet = None,
              key_outlet = None,
              step_callback = None,
              endcycle_callback = None,
              heat_transfer_callback = None,
              lump_energy_balance_callback = None,
              valves_callback = None,
              solver_method = 'Euler',
              OneCycle = False,
              Abort = None,
              pipe_abort = None,
              UseNR = True,
              **kwargs):
        """
        This is the driving function for the PDSim model.  It can be extended through the 
        use of the callback functions
        
        Parameters
        ----------
        key_inlet : string
            The key for the flow node that represents the upstream quasi-steady point
        key_outlet : string
            The key for the flow node that represents the upstream quasi-steady point
        step_callback : function
            This callback is passed on to the cycle() function. 
            See :func:`PDSim.core.Core.PDSimCore.cycle` for a description of the callback
        endcycle_callback : function
            This callback gets called at the end of the cycle to determine whether the cycle 
            should be run again.  It returns a flag(``redo`` that is ``True`` if the cycle 
            should be run again, or ``False`` if the cycle iterations have reached convergence.).  
            This callback does not take any inputs
        """
        
        # Set up a pipe for accepting a value of True which will abort the run
        # Typically used from the GUI
        self.pipe_abort = pipe_abort
        
        #If a function called pre_solve is provided, call it with no input arguments
        if hasattr(self,'pre_solve'):
            self.pre_solve()
            
        #This runs before the model starts at all
        self.__pre_run()
        
        #Both inlet and outlet keys must be connected to invariant nodes - 
        # that is they must be part of the tubes which are all quasi-steady
        if not key_inlet == None and not key_inlet in self.Tubes.Nodes:
            raise KeyError('key_inlet must be a Tube node')
        if not key_outlet == None and not key_outlet in self.Tubes.Nodes:
            raise KeyError('key_outlet must be a Tube node')
        
        self.key_inlet = key_inlet
        self.key_outlet = key_outlet
        
        if Abort is None and pipe_abort is not None:
            Abort = self.check_abort
        
        # (2) Run a cycle with the given values for the temperatures
        self.__pre_cycle()
        #x_state only includes the values for the chambers, the valves start closed
        #Since indexed, yields a copy
        self.x_state = self.__get_from_matrices(0)[0:len(self.stateVariables)*self.CVs.N]
            
        def OBJECTIVE_ENERGY_BALANCE(Td_Tlumps):
            print Td_Tlumps,'Td,Tlumps inputs'
            # Td_Tlumps is a list (or listm or np.ndarray)
            Td_Tlumps = list(Td_Tlumps)
            # Consume the first element as the discharge temp 
            self.Td = float(Td_Tlumps.pop(0))
            # The rest are the lumps in order
            self.Tlumps = Td_Tlumps
            
            # (0). Update the discharge temperature
            p = self.Tubes.Nodes[key_outlet].p
            self.Tubes.Nodes[key_outlet].update({'T':self.Td,'P':p})
            
            # (1). First, run all the tubes
            for tube in self.Tubes:
                tube.TubeFcn(tube)
            
            def OBJECTIVE_CYCLE(x_state):
                #Convert numpy array to listm
                x_state = listm(x_state)
                print 'x_state is', x_state
                
                t1 = clock()
                if solver_method == 'Euler':
                    if hasattr(self,'EulerN'):
                        N=self.EulerN
                    else:
                        N=7000
                    self.cycle_SimpleEuler(N,x_state,step_callback=step_callback,
                                           heat_transfer_callback=heat_transfer_callback,
                                           valves_callback=valves_callback)
                elif solver_method == 'Heun':
                    if hasattr(self,'HeunN'):
                        N=self.HeunN
                    else:
                        N=7000
                    self.cycle_SimpleEuler(N,x_state,step_callback=step_callback,
                                           heat_transfer_callback=heat_transfer_callback,
                                           valves_callback=valves_callback)
                elif solver_method == 'RK45':
                    if hasattr(self,'RK45_eps'):
                        eps_allowed=self.RK45_eps
                    else:
                        eps_allowed = 1e-8
                    self.cycle_RK45(x_state,
                                    eps_allowed = eps_allowed,
                                    step_callback=step_callback,
                                    heat_transfer_callback=heat_transfer_callback,
                                    valves_callback=valves_callback,
                                    **kwargs)
                else:
                    raise AttributeError
                t2 = clock()
                print 'Elapsed time for cycle is {0:g} s'.format(t2-t1)
                
                if (key_inlet in self.FlowsProcessed.mean_mdot and 
                    key_outlet in self.FlowsProcessed.mean_mdot):
                    
                    mdot_out = self.FlowsProcessed.mean_mdot[key_outlet]
                    mdot_in = self.FlowsProcessed.mean_mdot[key_inlet]
                    print 'Mass flow difference',(mdot_out+mdot_in)/mdot_out*100,' %'
                    
                    h_outlet = (self.FlowsProcessed.integrated_mdoth[key_outlet]
                                /self.FlowsProcessed.integrated_mdot[key_outlet])
                    p = self.Tubes.Nodes[key_outlet].p
                    Fluid = self.Tubes.Nodes[key_outlet].Fluid
                    # Multimply the residual on the discharge temperature by the capacitance rate
                    # to get all the residuals in units of kW
                    C_outlet = self.Tubes.Nodes[key_outlet].cp*mdot_out
                    self.resid_Td = C_outlet*(Props('T','H',h_outlet,'P',p,Fluid) - self.Td)
                else:
                    raise KeyError
                      
                #debug_plots(self)
                
                #If the abort function returns true, quit this loop
                if (Abort is not None and Abort()) or OneCycle:
                    print 'Quitting OBJECTIVE_CYCLE loop in core.solve'
                    return None #Stop
                        
                if endcycle_callback is None:
                    return None #Stop
                else:
                    #endcycle_callback returns the errors and new initial st
                    errors, x_state_ = endcycle_callback()
                    self.x_state = x_state_[:] #Make a copy
                    return errors
                
            #! End of OBJECTIVE_CYCLE
            
            if UseNR:
                self.x_state = Broyden(OBJECTIVE_CYCLE, self.x_state, Nfd = 1, dx = 0.01*np.array(self.x_state), itermax = 50, ytol = 1e-4)
                print self.x_state
                if self.x_state[0] is None:
                    return None
            else:
                diff_old = None
                x_state_prior = None
                
                while True:
                    x_state_old = listm(self.x_state).copy()
                    OBJECTIVE_CYCLE(x_state_old)
                    x_state_new = listm(self.x_state).copy()
                    
                    try:    
                        diff = (x_state_old - x_state_new) / x_state_old
                        diff_abs = [abs(dx) for dx in (x_state_old - x_state_new)/x_state_old]
                    except IndexError:
                        continue
                    
                    def all_negative(yvec):
                        for y in yvec:
                            if y>0:
                                return False
                        return True
                    
                    print 'diff', diff
                    if max(diff_abs) < 1e-4:
                        break
#                    elif diff_old is not None and all_negative(diff * diff_old):
#                        #If the signs are opposite for all of the diffs, weight the next step
#                        print 'Weighting old and new initial states'
#                        print 'new', diff, x_state_old
#                        print 'old', diff_old, x_state_prior
#                        self.x_state = (x_state_prior + x_state_old)/2
                    else:
                        self.x_state = x_state_new
                    
#                    diff_old = diff.copy()
#                    x_state_prior = x_state_old.copy()
                        
                    
            # (3) After convergence of the inner loop, check the energy balance on the lumps
            if lump_energy_balance_callback is not None:
                resid_HT = lump_energy_balance_callback()
                if not isinstance(resid_HT,list) and not isinstance(resid_HT,listm):
                    resid_HT = [resid_HT]
                
            #If the abort function returns true, quit this loop
            if (Abort is not None and Abort()) or OneCycle:
                print 'Quitting OBJECTIVE function in core.solve'
                return None
                
            print [self.resid_Td]+resid_HT,'resids'
            return [self.resid_Td]+resid_HT
        
        #end of OJECTIVE
        
        x_soln = Broyden(OBJECTIVE_ENERGY_BALANCE,[360.0,390.0],dx=0.3,ytol=0.001,itermax=30)
        print 'Solution is', x_soln,'Td, Tlumps'
        
        del self.resid_Td, self.x_state
        if Abort is None or not Abort():
            self.__post_solve()
        
    def __post_solve(self):
        """
        Do some post-processing to calculate flow rates, efficiencies, etc.  
        """
        if not hasattr(self,'Qamb'):
            self.Qamb = 0
            
        #The total mass flow rate
        self.mdot = self.FlowsProcessed.mean_mdot[self.key_inlet]
        
        for key, State in self.Tubes.Nodes.iteritems():
            if key == self.key_inlet:
                 inletState = State
            if key == self.key_outlet:
                 outletState = State
        
        if callable(self.Vdisp):
            Vdisp = self.Vdisp()
        else:
            Vdisp = self.Vdisp
            
        self.eta_v = self.mdot / (self.omega/(2*pi)*Vdisp*inletState.rho)
        h1 = inletState.h
        h2 = outletState.h
        s1 = inletState.s
        print h1, h2
        T2s = newton(lambda T: Props('S','T',T,'P',outletState.p,outletState.Fluid)-s1,inletState.T+30)
        h2s = Props('H','T',T2s,'P',outletState.p,outletState.Fluid)
        self.eta_a = (h2s-h1)/(h2-h1)
        self.Wdot_i = self.mdot*(h2s-h1)
        # self.Qamb is positive if heat is being added to the lumped mass
        self.Wdot = self.mdot*(h2-h1)-self.Qamb
        
        #Resize all the matrices to keep only the real data
        self.t = self.t[  0:self.Ntheta]
        self.T = self.T[:,0:self.Ntheta]
        self.p = self.p[:,0:self.Ntheta]
        self.m = self.m[:,0:self.Ntheta]
        self.rho = self.rho[:,0:self.Ntheta]
        self.V = self.V[:,0:self.Ntheta]
        self.dV = self.dV[:,0:self.Ntheta]
        self.xValves = self.xValves[:,0:self.Ntheta]
        
        print 'mdot*(h2-h1),P-v,Qamb', self.Wdot, self.Wdot_pv, self.Qamb
        print 'Mass flow rate is',self.mdot*1000,'g/s'
        print 'Volumetric efficiency is',self.eta_v*100,'%'
        print 'Adiabatic efficiency is',self.eta_a*100,'%'
        
    def derivs(self,theta,x,heat_transfer_callback=None,valves_callback=None):
        """
        
        derivs() is an internal function that should (probably) not actually be called by
        any user-provided code, but is documented here for completeness.
        
        Parameters
        ----------
        theta : float
            The value of the independent variable
        x : ``listm`` instance
            The array of the state variables (plus valve parameters) 
        heat_transfer_callback : ``None`` or function
            A function that has the form::
            
                Q=heat_transfer_callback(theta)
            
            where ``Q`` is the ``listm`` instance of heat transfer rates to the control 
            volume.  If an entry of ``Q`` is positive, heat is being transferred 
            TO the fluid contained in the control volume
            
            If ``None``, no heat transfer is used
                
        Returns
        -------
        dfdt : ``listm`` instance
        
        """
        
        #Call the Cython method
        #return self._derivs(theta,x,heat_transfer_callback)
        
        #1. Calculate the volume and derivative of volume of each control volumes
        #    Return two lists, one for the volumes, second for the derivatives of volumes 
        V,dV=self.CVs.volumes(theta)

        if self.__hasLiquid__==True:
            raise NotImplementedError
        else:
            self.CVs.updateStates('T',x[0:self.CVs.Nexist],'D',x[self.CVs.Nexist:2*self.CVs.Nexist])
        
        ## Calculate properties and property derivatives
        ## needed for differential equations
        rho = listm(self.CVs.rho)
        v = 1.0/rho
        m = V*rho
        T = listm(self.CVs.T)
        h = listm(self.CVs.h)
        cv = listm(self.CVs.cv)
        dpdT = listm(self.CVs.dpdT)
        p = listm(self.CVs.p)
        
        self.V_=V
        self.dV_=dV
        self.h_ = h
        self.rho_=rho
        self.m_=m
        self.p_=p
        self.T_=T

        #Run the valves model if it is provided in order to calculate the new valve positions
        if self.__hasValves__==True:
            for iV,Valve in enumerate(self.Valves):
                Valve.set_xv(x[2*self.CVs.Nexist+iV*2:2*self.CVs.Nexist+(iV+1)*2])
            f_valves=valves_callback()
        
        #Build a dictionary of enthalpy values to speed up _calculate
        #Halves the number of calls to get_h which calculates the enthalpy
        hdict = {key:h_ for key, h_ in zip(self.CVs.exists_keys, h)}
        hdict.update(self.Tubes_hdict)
        
        #2. Calculate the mass flow terms between the model components
        self.Flows.calculate(self, hdict)
        summerdT,summerdm=self.Flows.sumterms(self)
        
        #3. Calculate the heat transfer terms
        if heat_transfer_callback is not None:
            Q=heat_transfer_callback(theta)
            if not len(Q) == self.CVs.Nexist:
                raise ValueError('Length of Q is not equal to length of number of CV')
        else:
            Q=0.0
        
        self.Q_=Q        
        dudxL=0.0
        summerdxL=0.0
        xL=0.0
        
        dmdtheta=summerdm
        dxLdtheta=1.0/m*(summerdxL-xL*dmdtheta)
        dTdtheta=1/(m*cv)*(-1.0*T*dpdT*(dV-v*dmdtheta)-m*dudxL*dxLdtheta-h*dmdtheta+Q/self.omega+summerdT)
        drhodtheta = 1.0/V*(dmdtheta-rho*dV)
        
        if self.__hasLiquid__==True:
            f=np.zeros((3*self.NCV,))
            f[0:self.NCV]=dTdtheta
            f[self.NCV:2*self.NCV]=dmdtheta
            f[2*self.NCV:3*self.NCV]=dxLdtheta
            return f
        else:
            f=list(dTdtheta)+list(drhodtheta)
            if self.__hasValves__==True:
                f+=f_valves
            return listm(f)
    
    def valves_callback(self):
        """
        This is the default valves_callback function that builds the list of 
        derivatives of position and velocity with respect to the crank angle
        
        It returns a ``list`` instance with the valve return values in order 
        """
        #Run each valve model in turn to calculate the derivatives of the valve parameters
        # for each valve
        f=[]
        for Valve in self.Valves:
            f+=Valve.derivs(self)
        return f
        
    def endcycle_callback(self,eps_wrap_allowed=0.0001):
        """
        This function can be called at the end of the cycle if so desired.
        Its primary use is to determine whether the cycle has converged for a 
        given set of discharge temperatures and lump temperatures.
        
        Parameters
        ----------
        eps_wrap_allowed : float
            Maximum error allowed, in absolute value
        
        Returns 
        -------
        redo : boolean
            ``True`` if cycle should be run again with updated inputs, ``False`` otherwise.
            A return value of ``True`` means that convergence of the cycle has been achieved
        """
        
        #debug_plots(self)
        
        #old and new CV keys
        LHS,RHS=[],[]
        errorT,error_rho,newT,new_rho,oldT,old_rho={},{},{},{},{},{}
        
        for key in self.CVs.exists_keys:
            # Get the 'becomes' field.  If a list, parse each fork of list. If a single key convert 
            # into a list so you can use the same code  
            if not isinstance(self.CVs[key].becomes,list):
                becomes=[self.CVs[key].becomes]
            else:
                becomes=self.CVs[key].becomes
            Iold = self.CVs.index(key)

            for newkey in becomes:
                Inew = self.CVs.index(newkey)
                newCV = self.CVs[newkey]
                # There can't be any overlap between keys
                if newkey in newT:
                    raise KeyError
                #What the temperature and density were at the start of the rotation
                oldT[newkey]=self.T[Inew, 0]
                old_rho[newkey]=self.rho[Inew, 0]
                #What they are at the end of the rotation
                newT[newkey]=self.T[Iold,self.Itheta-1]
                new_rho[newkey]=self.rho[Iold,self.Itheta-1]
                errorT[newkey]=(oldT[newkey]-newT[newkey])/newT[newkey]
                error_rho[newkey]=(old_rho[newkey]-new_rho[newkey])/new_rho[newkey]
                #Update the list of keys for setting the exist flags
                LHS.append(key)
                RHS.append(newkey)
                #Update the CV with the new value - use the copy of the CVs from before
                #newCV.State.update({'T':self.T[Iold,self.Itheta],'D':self.rho[Iold,self.Itheta]})
            
        errorTlist = [errorT[key] for key in self.CVs.keys() if key in newT]
        errorrholist = [error_rho[key] for key in self.CVs.keys() if key in new_rho]
        newTlist = [newT[key] for key in self.CVs.keys() if key in newT]
        new_rho_list = [new_rho[key] for key in self.CVs.keys() if key in new_rho]
        errorT_abs = [abs(err) for err in errorTlist]
        error_rho_abs = [abs(err) for err in errorrholist]
        
        #Reset the exist flags for the CV - this should handle all the possibilities
        #Turn off the LHS CV
        for key in LHS:
            self.CVs[key].exists=False
        #Turn on the RHS CV
        for key in RHS:
            self.CVs[key].exists=True
            
        self.CVs.rebuild_exists()
        
        return errorTlist + errorrholist, newTlist + new_rho_list
    
if __name__=='__main__':
    print 'This is the base class that is inherited by other compressor types.  Running this file doesn\'t do anything'
