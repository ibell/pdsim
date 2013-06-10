from __future__ import division

from math import pi
import textwrap
import copy
from time import clock
from cPickle import loads, dumps
import inspect

##--  Package imports  --
from PDSim.flow import flow,flow_models
from containers import STATE_VARS_TM, CVArrays
from PDSim.flow.flow import FlowPathCollection, FlowPath
from containers import ControlVolumeCollection,Tube,TubeCollection
from PDSim.plot.plots import debug_plots
from PDSim.misc.solvers import Broyden, MultiDimNewtRaph
from PDSim.misc.datatypes import arraym, empty_arraym
from _core import delimit_vector, setcol, getcol, _PDSimCore
import PDSim.core.callbacks

##-- Non-package imports  --
import numpy as np
from scipy.integrate import trapz, simps
from scipy.optimize import newton
import h5py

## matplotlib imports
import matplotlib.pyplot as plt
import pylab

## Coolprop imports
from CoolProp.CoolProp import Props
from CoolProp.State import State

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
        return func.__get__(obj, cls)

import copy_reg
import types
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

#An empty class for storage
class struct(object):
    pass    

class PDSimCore(_PDSimCore):
    """
    This is the main driver class for the model
    
    This class is not intended to be run on its own.  It must be subclassed and extended to provide functions for mass flow, etc. 
    
    The basic order of steps that should be followed can be summarized as
    #. Instantiate the subclass of PDSimCore
    #. Add each of the control volumes
    #. Add each of the tubes
    #. Add all the flow models between CV and tubes
    #. Add valves (if applicable)
    #. Connect the callbacks for heat transfer, step, etc.
    #. Run the model
    
    """
    def __init__(self,stateVariables=None):
        """
        Initialization of the PDSimCore
        
        Parameters
        ----------
        stateVariables : mutable object [list or tuple], optional
            list of keys for the state variables to be used.  Current options are 'T','D' or 'T','M'.  Default state variables are 'T','M'
        """
        #Initialize the containers to be empty
        
        #: The Valves container class
        self.Valves=[]
        #: The :class:`ControlVolumeCollection <PDSim.core.containers.ControlVolumeCollection>` instance
        #: that contains all the control volumes in the machine
        self.CVs=ControlVolumeCollection()
        #: The :class:`FlowPathCollection <PDSim.flow.flow.FlowPathCollection>` 
        #: instance
        self.Flows=FlowPathCollection()
        #: A :class:`list` that contains copies of the 
        #: :class:`FlowPathCollection <PDSim.flow.flow.FlowPathCollection>` 
        #: at each crank angle
        self.FlowStorage=[]
        
        self.Tubes=TubeCollection()
        self.Tlumps=np.zeros((1,1))
        self.steps=[]
        self.__hasValves__=False
        
        #: A storage of the initial state vector
        self.xstate_init = None
        
        #: A storage of the initial valves vector
        if isinstance(stateVariables,(list,tuple)):
            self.stateVariables=list(stateVariables)
        else:
            self.stateVariables=['T','M']
        self._want_abort = False
        
        #Build a structure to hold all the callbacks
        self.callbacks = PDSim.core.callbacks.CallbackContainer()
    
    def _get_from_matrices(self,i):
        """
        Get values back from the matrices and reconstruct the state variable list
        """
        
        if self.__hasLiquid__==True:
            raise NotImplementedError
            #return np.hstack([self.T[:,i],self.m[:,i],self.xL[:,i]])
        else:
            VarList=[]
            exists_indices = self.CVs.exists_indices
            for s in self.stateVariables:
                if s=='T':
                    VarList+=getcol(self.T,i,exists_indices)
                elif s=='D':
                    VarList+=getcol(self.rho,i,exists_indices)
                elif s=='M':
                    VarList+=getcol(self.m,i,exists_indices)
                else:
                    raise KeyError
            if self.__hasValves__:
                VarList+=list(self.xValves[:,i])
                
            return arraym(VarList)
        
    def _statevars_to_dict(self,x):
        d={}
        for iS,s in enumerate(self.stateVariables):
            x_=(x[iS*self.CVs.Nexist:self.CVs.Nexist*(iS+1)])
            if s=='T':
                d['T']=x_
            elif s=='D':
                d['D']=x_
            elif s=='M':
                d['M']=x_
        return d
        
    def _put_to_matrices(self,x,i):
        """
        Take a state variable list and put back in numpy matrices
        """
        exists_indices=self.CVs.exists_indices
        Nexist = self.CVs.Nexist
        Ns = len(self.stateVariables)
        if self.__hasLiquid__==True:
            raise NotImplementedError
            self.T[:,i]=x[0:self.NCV]
            self.m[:,i]=x[self.NCV:2*self.NCV]
            self.xL[:,i]=x[2*self.NCV:3*self.NCV]
        else: # self.__hasLiquid__==False
            for iS, s in enumerate(self.stateVariables):
                x_=arraym(x[iS*self.CVs.Nexist:self.CVs.Nexist*(iS+1)])
                if s=='T':
                    setcol(self.T, i, exists_indices, x_)
                elif s=='D':
                    setcol(self.rho, i, exists_indices, x_)
                elif s=='M':
                    setcol(self.m, i, exists_indices, x_)
            #Left over terms are for the valves
            if self.__hasValves__:
                setcol(self.xValves, i, 
                       range(len(self.Valves)*2), 
                       arraym(x[Ns*Nexist:len(x)])
                       )
            
            # In the first iteration, self.core has not been filled, so do not 
            # overwrite with the values in self.core.m and self.core.rho
            if self.core.m[0] > 0.0 :
                setcol(self.m, i, exists_indices, self.core.m)
            if self.core.rho[0] > 0.0 :
                setcol(self.rho, i, exists_indices, self.core.rho)
            
            setcol(self.V, i, exists_indices, self.core.V)
            setcol(self.dV, i, exists_indices, self.core.dV)
            setcol(self.p, i, exists_indices, self.core.p)
            setcol(self.h, i, exists_indices, self.core.h)    
            setcol(self.Q, i, exists_indices, self.core.Q)
    
    def _postprocess_flows(self):
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
            return flow.sumterms_given_CV(key, Flows)
            
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
            self.FlowsProcessed.summed_mdot[key]=np.zeros((self.Ntheta,))
            self.FlowsProcessed.summed_mdoth[key]=np.zeros((self.Ntheta,))
            
            assert self.Ntheta == len(self.FlowStorage)
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
        # The other one is invariant because it is quasi-steady.
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
            
            #For each tube, update the flow going through it
            #Tube.mdot is always a positive value
            Tube.mdot = max(abs(mdot1), abs(mdot2))
            
        self.mdot = self.FlowsProcessed.mean_mdot[self.key_inlet]
            
    def _postprocess_HT(self):
        """
        Postprocess the heat transfer terms
        
        Here we
        calculate the mean heat transfer rate over the course of the cycle 
        """
        self.HTProcessed=struct()
        r = range(self.Ntheta)
        
        #Remove all the NAN placeholders and replace them with zero values
        self.Q[np.isnan(self.Q)] = 0.0
        
        #Sum at each step of the revolution
        self.HTProcessed.summed_Q = np.sum(self.Q, axis = 0) #kW
        
        #Get the mean heat transfer rate
        self.HTProcessed.mean_Q = trapz(self.HTProcessed.summed_Q[r], self.t[r])/(self.t[self.Ntheta-1]-self.t[0])
    
    def guess_outlet_temp(self, inletState, p_outlet, eta_a = 0.7):
        """ 
        Function to guess outlet temperature
        
        Using a guess value for the adiabatic efficiency, calculate the guessed
        outlet temperature.  In compressor mode, the adiabatic efficiency is defined by
        
        .. math::
        
            \eta_a = \\frac{h_{2s}-h_1}{h_2-h_1}
            
        and in expander mode it is defined by
        
        .. math::
        
            \eta_a = \\frac{h_2-h_1}{h_{2s}-h_1}
            
        This function can also be overloaded by the subclass in order to 
        implement a different guess method
        """
        
        h1 = inletState.h
        h2s = Props('H','S',inletState.s,'P',p_outlet, inletState.Fluid)
        if p_outlet > inletState.p:
            #Compressor Mode
            h2 = h1 + (h2s-h1)/eta_a
            return Props('T','H',h2,'P',p_outlet, inletState.Fluid)
        else:
            #Expander Mode
            h2 = h1 + (h2s-h1)*eta_a
            return Props('T','H',h2,'P',p_outlet, inletState.Fluid)
    
    def reset_initial_state(self):
        """
        Reset the initial state of the core class, typically after doing a 
        preconditioning run
        """
        
        for k,CV in zip(self.CVs.keys,self.CVs.CVs):
            if k in self.exists_CV_init:
                CV.exists = True
            else:
                CV.exists = False

        #Update the existence of each of the CV
        self.update_existence()
             
        #Only the State variables, not the valves
        self.x_state = self.xstate_init
        #Make a copy
        x = self.xstate_init.copy()
        #Add the values from the valves
        if self.__hasValves__:
            x.extend(empty_arraym(2*len(self.Valves)))
        self._put_to_matrices(x, 0)
        #Reset the temporary variables
        self.xstate_init = None
        self.exists_CV_init = None
                
    def update_existence(self):
        """
        Update existence flags for Tubes and control volumes
        
        This function is required to be called when the existence of any control
        volume or tube changes in order to ensure that internal flags are set
        properly
        """
        
        # Update the existence flags in all the control volumes
        self.CVs.rebuild_exists()
        
        # Update the array of enthalpies in the tubes
        self.Tubes.update_existence(self.CVs.Nexist)
        
        # Update the existence of each of the flows
        self.Flows.update_existence(self)
        
        # Update the sizes of the internal arrays in self.core  
        self.core.update_size(self.CVs.Nexist)
        
    def add_flow(self,FlowPath):
        """
        Add a flow path to the model
        
        Parameters
        ----------
        FlowPath : :class:`FlowPath <PDSim.flow.flow.FlowPath>`  instance
            An initialized flow path 
        """
        #Add FlowPath instance to the list of flow paths
        self.Flows.append(FlowPath)
        
    def add_CV(self,CV):
        """
        Add a control volume to the model
        
        Parameters
        ----------
        CV : :class:`ControlVolume <PDSim.core.containers.ControlVolume>` instance
            An initialized control volume
        """
        
        if CV.key in self.CVs.keys:
            raise KeyError('Sorry but the key for your Control Volume ['+CV.key+'] is already in use')
        
        #Add the CV to the collection
        self.CVs.add(CV)
        self.CVs.rebuild_exists()
        
    def add_tube(self,Tube):
        """
        Add a tube to the model.
        
        Parameters
        ----------
        Tube : :class:`Tube <PDSim.core.containers.Tube>` instance
            An initialized tube.
        """
        #Add it to the list
        self.Tubes.append(Tube)
        self.Tubes.update()
        
    def add_valve(self,Valve):
        """
        Add a valve to the model.
        
        Parameters
        ----------
        Valve : :class:`ValveModel <PDSim.flow.flow_models.ValveModel>` instance
            An initialized valve.
        """
        #Add it to the list
        self.Valves.append(Valve)
        self.__hasValves__=True
        
    def pre_run(self, N = 40000):
        """
        This function gets called before the run begins.  It builds large matrices
        to store values, and does other initialization. 
        """
        # Build the full numpy arrays for temperature, volume, etc.
        self.t=np.zeros((N,))
        self.T=np.zeros((self.CVs.N,N))
        self.T.fill(np.nan)
        self.p=self.T.copy()
        self.h = self.T.copy()
        self.m = self.T.copy()
        self.V = self.T.copy()
        self.dV = self.T.copy()
        self.rho = self.T.copy()
        self.Q = self.T.copy()
        self.xL = self.T.copy()
        self.xValves = np.zeros((2*len(self.Valves),N))
        
        # Initialize the core class that contains the arrays and the derivs
        self.core = CVArrays(0)
        
        # Update the existence of all the control volumes
        self.update_existence()
        
        # Set a flag about liquid flooding
        self.__hasLiquid__ = False
        
    def pre_cycle(self, x0 = None):
        """
        This runs before the cycle is run but after pre_run has been called
        
        Parameters
        ----------
        x0 : :class:`arraym <PDSim.misc.datatypes.arraym>` instance
        """        
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
        VdV=[CV.V_dV(0.0,**CV.V_dV_kwargs) for CV in self.CVs.exists_CV]
        V,dV = zip(*VdV)
        
        self.t[0]=0
        
        # If x0 is provided, use its values to initialize the chamber states
        if x0 is None:
            # self.CVs.exists_indices is a list of indices of the CV with the same order of entries
            # as the entries in self.CVs.T
            self.T[self.CVs.exists_indices, 0] = self.CVs.T
            self.p[self.CVs.exists_indices, 0] = self.CVs.p
            self.rho[self.CVs.exists_indices, 0] = self.CVs.rho
            self.m[self.CVs.exists_indices, 0] = self.CVs.rho*arraym(V)
        else:
            #x0 is provided, but need to pad it out to include valve values
            x0_ = x0.copy()
            
            # If x0 is longer than the product of the number of state variables 
            # and CV in existence, the valve data is already included and must not be 
            # added to the array of independent variables
            if self.__hasValves__ and len(x0) == self.CVs.Nexist*len(self.stateVariables):
                #Load up the rest of the array with zeros since the valves start closed and at rest
                x0_.extend(empty_arraym(len(self.Valves)*2))
            self._put_to_matrices(x0_, 0)
        
        # Assume all the valves to be fully closed and stationary at the beginning of cycle
        self.xValves[:,0]=0
        
        self.Tubes_hdict={}
        for Tube in self.Tubes:
            self.Tubes_hdict[Tube.key1]=Tube.State1.get_h()
            self.Tubes_hdict[Tube.key2]=Tube.State2.get_h()
        
    def cycle_SimpleEuler(self,N,x_state,tmin=0,tmax=2*pi):
        """
        The simple Euler PDSim ODE integrator
        
        Parameters
        ----------
        N : integer
            Number of steps taken.  There will be N+1 entries in the state matrices
        x_state : 
            The initial values of the variables (ONLY the state variables, no valves)
        tmin : float, optional
            Starting value of the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``]
        tmax : float, optional
            Ending value for the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``] 
        
        """
        #Do some initialization - create arrays, etc.
        self.pre_cycle(x_state)
        
        #Step variables
        t0=tmin
        h=(tmax-tmin)/(N)
        
        # Get the beginning of the cycle configured
        # Put a copy of the values into the matrices
        xold = x_state.copy()
        #Add zeros for the valves as they are assumed to start closed and at rest
        if self.__hasValves__:
            xold.extend(empty_arraym(len(self.Valves)*2))
        self._put_to_matrices(xold, 0)
        
        for Itheta in range(N):
            #Once every 100 steps check if you are supposed to abort
            if self._check_cycle_abort(Itheta):
                return 'abort'
                  
            #Call the step callback if provided
            if self.callbacks.step_callback is not None:
                h = self.callbacks.step_callback(t0, h, Itheta)
                disable = self.callbacks.step_callback.disable_adaptive
                if disable:
                    print 'CV have changed'
                    xold = self._get_from_matrices(Itheta)
            
            # Step 1: derivatives evaluated at old values of t = t0
            f1 = self.derivs(t0, xold)
            xnew = xold+h*f1
            
            #Store at the current index (first run at index 0)
            self.t[Itheta] = t0
            self._put_to_matrices(xold, Itheta)
            self.FlowStorage.append(self.Flows.get_deepcopy())
            
            # Everything from this step is finished, now update for the next
            # step coming
            t0+=h
            xold = xnew
            
        #Run this at the end at index N-1
        #Run this at the end
        V,dV=self.CVs.volumes(t0)
        #Stored at the old value
        self.t[N]=t0
        
        self._put_to_matrices(xnew, N)
        
        #ensure you end up at the right place
        assert abs(t0-tmax)<1e-10
        
        self.derivs(t0,xold)
        self.FlowStorage.append(self.Flows.get_deepcopy())
        
        if sorted(self.stateVariables) == ['D','T']:
            self.CVs.updateStates('T',xnew[0:self.CVs.Nexist],'D',xnew[self.CVs.Nexist:2*self.CVs.Nexist])
        elif sorted(self.stateVariables) == ['M','T']:
            self.CVs.updateStates('T',xnew[0:self.CVs.Nexist],'D',xnew[self.CVs.Nexist:2*self.CVs.Nexist]/V)
        else:
            raise NotImplementedError
        
        # last index is Itheta, number of entries in FlowStorage is Ntheta
        print 'Number of steps taken', N
        self.Itheta = N
        self.Ntheta = N+1
        self.post_cycle()
        
    def cycle_Heun(self, N, x_state, tmin = 0, tmax = 2*pi):
        """
        Use the Heun method (modified Euler method)
        
        Parameters
        ----------
        N : integer
            Number of steps to take (N+1 entries in the state vars matrices)
        x_state : 
            The initial values of the variables (only the state variables)
        tmin : float
            Starting value of the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``]
        tmax : float
            Ending value for the independent variable.  ``t`` is in the closed range [``tmin``, ``tmax``] 
        
        """
        #Do some initialization - create arrays, etc.
        self.pre_cycle()
        
        #Start at an index of 0
        Itheta=0
        t0=tmin
        h=(tmax-tmin)/(N)
        
        # Get the beginning of the cycle configured
        # Put a copy of the values into the matrices
        self._put_to_matrices(x_state.copy(), 0)
    
        for Itheta in range(N):
            
            #Once every 100 steps check if you are supposed to abort
            if self._check_cycle_abort(Itheta):
                return 'abort'
            
            if self.callbacks.step_callback!=None:
                self.callbacks.step_callback(t0,h,Itheta)
                
            xold=self._get_from_matrices(Itheta)
                        
            # Step 1: derivatives evaluated at old values
            f1=self.derivs(t0,xold)
            xtemp=xold+h*f1
            
            #Stored at the starting value of the step
            self.t[Itheta]=t0+h
            self._put_to_matrices(xold,Itheta)
            self.FlowStorage.append(self.Flows.get_deepcopy())
            
            #Step 2: Evaluated at predicted step
            f2=self.derivs(t0+h,xtemp)
            xnew = xold + h/2.0*(f1 + f2)
            
            t0+=h
            xold = xnew
            
        #ensure you end up at the right place
        assert abs(t0-tmax)<1e-10
        
        #Run this at the end
        V,dV=self.CVs.volumes(t0)
        #Stored at the old value
        self.t[N]=t0
        self.derivs(t0,xold)
        self._put_to_matrices(xnew,N)
        self.FlowStorage.append(self.Flows.get_deepcopy())
        if sorted(self.stateVariables) == ['D','T']:
            self.CVs.updateStates('T',xnew[0:self.CVs.Nexist],'D',xnew[self.CVs.Nexist:2*self.CVs.Nexist])
        elif sorted(self.stateVariables) == ['M','T']:
            self.CVs.updateStates('T',xnew[0:self.CVs.Nexist],'D',xnew[self.CVs.Nexist:2*self.CVs.Nexist]/V)
        else:
            raise NotImplementedError
        
        print 'Number of steps taken', N,'len(FlowStorage)',len(self.FlowStorage)
        self.Itheta = N
        self.Ntheta = N+1
        self.post_cycle()
        return
        
    def cycle_RK45(self,x_state,tmin=0,tmax=2.0*pi,hmin=1e-4,eps_allowed=1e-10,step_relax=0.9,
              valves_callback = None,
              UseCashKarp=True,**kwargs):
        """
        
        This function implements an adaptive Runge-Kutta-Feldberg 4th/5th order
        solver for the system of equations
        
        Parameters
        ----------
        x_state : arraym
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
        self.pre_cycle()
        
        #Start at an index of 0
        Itheta = 0
        t0 = tmin
        h = hmin
        
        # Get the beginning of the cycle configured
        # Put a copy of the values into the matrices
        xold = x_state.copy()
        #Add zeros for the valves as they are assumed to start closed and at rest
        if self.__hasValves__:
            xold.extend(empty_arraym(len(self.Valves)*2))
        self._put_to_matrices(xold, 0)
        
        gamma1=16.0/135.0
        gamma2=0.0
        gamma3=6656.0/12825.0
        gamma4=28561.0/56430.0
        gamma5=-9.0/50.0
        gamma6=2.0/55.0
        
        #t is the independent variable here, where t takes on values in the bounded range [tmin,tmax]
        while (t0<tmax):
            
            #Once every 100 steps check if you are supposed to abort
            if self._check_cycle_abort(Itheta):
                return 'abort'
            
            stepAccepted=False
            while not stepAccepted:
                
                #Reset the flag
                disableAdaptive=False
                
                if t0 + h > tmax:
                    disableAdaptive = True
                    h = tmax - t0
            
                if self.callbacks.step_callback is not None and not disableAdaptive:
                    #The user has provided a disabling function for the adaptive method
                    #Call it to figure out whether to use the adaptive method or not
                    #Pass it a copy of the compressor class and the current step size
                    #The function can modify anything in the class to change flags, existence, merge, etc.
                    h = self.callbacks.step_callback(t0, h, Itheta)
                    disableAdaptive = self.callbacks.step_callback.disable_adaptive
                    x = self._get_from_matrices(Itheta)
                    
                    if disableAdaptive and np.all(np.isfinite(x)):
                        xold = self._get_from_matrices(Itheta)
                else:
                    disableAdaptive=False
                    
                if h < hmin and not disableAdaptive:
                    #Step is too small, just use the minimum step size
                    h = 1.0*hmin
                    disableAdaptive = True
                    
                if not UseCashKarp:
                    ## Using RKF constants ##
                    
                    # Step 1: derivatives evaluated at old values
                    f1=self.derivs(t0,xold)
                    xnew1=xold+h*(+1.0/4.0*f1)
                    
                    #Store a copy of the flows for future use as well as a buffered set of state variables
                    Flows_temporary = self.Flows.get_deepcopy()
                    core_temporary = self.core.copy()
                    
                    f2=self.derivs(t0+1.0/4.0*h,xnew1)
                    xnew2=xold+h*(+3.0/32.0*f1+9.0/32.0*f2)
    
                    f3=self.derivs(t0+3.0/8.0*h,xnew2)
                    xnew3=xold+h*(+1932.0/2197.0*f1-7200.0/2197.0*f2+7296.0/2197.0*f3)
    
                    f4=self.derivs(t0+12.0/13.0*h,xnew3)
                    xnew4=xold+h*(+439.0/216.0*f1-8.0*f2+3680.0/513.0*f3-845.0/4104.0*f4)
                    
                    f5=self.derivs(t0+h,xnew4)
                    xnew5=xold+h*(-8.0/27.0*f1+2.0*f2-3544.0/2565.0*f3+1859.0/4104.0*f4-11.0/40.0*f5)
                    
                    #Updated values at the next step
                    f6=self.derivs(t0+h/2.0,xnew5)
                    
                    xnew=xold+h*(gamma1*f1 + gamma2*f2 + gamma3*f3 + gamma4*f4 + gamma5*f5 + gamma6*f6)
                        
                    error=h*(1.0/360.0*f1-128.0/4275.0*f3-2197.0/75240.0*f4+1.0/50.0*f5+2.0/55.0*f6)
                else:
                    # Step 1: derivatives evaluated at old values
                    f1=self.derivs(t0,xold)
                    xnew1=xold+h*(+1.0/5.0*f1)
                    
                    #Store a copy of the flows for future use as well as a buffered set of state variables
                    Flows_temporary = self.Flows.get_deepcopy()
                    core_temporary = self.core.copy()
                    
                    f2=self.derivs(t0+1.0/5.0*h,xnew1)
                    xnew2=xold+h*(+3.0/40.0*f1+9.0/40.0*f2)
    
                    f3=self.derivs(t0+3.0/10.0*h,xnew2)
                    xnew3=xold+h*(3.0/10.0*f1-9.0/10.0*f2+6.0/5.0*f3)
    
                    f4=self.derivs(t0+3.0/5.0*h,xnew3)
                    xnew4=xold+h*(-11.0/54.0*f1+5.0/2.0*f2-70/27.0*f3+35.0/27.0*f4)
                    
                    f5=self.derivs(t0+h,xnew4)
                    xnew5=xold+h*(1631.0/55296*f1+175.0/512.0*f2+575.0/13824.0*f3+44275.0/110592.0*f4+253.0/4096.0*f5)
                    
                    f6=self.derivs(t0+7/8*h,xnew5)
                    
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
                        h *= step_relax*(eps_allowed/max_error)**(0.3)
#                        print h, eps_allowed, max_error
                        stepAccepted=False
                    else:
                        # Accept the step regardless of whether the error 
                        # is too large or not
                        stepAccepted = True
                else:
                    stepAccepted = True
                
            #This block is for saving of values
            #
            #Store crank angle at the current index (first run at Itheta=0)
            self.t[Itheta] = t0

            # Use the copy that was stored before
            self.core = core_temporary.copy()
            
            # Store the values for volumes and state vars in the matrices
            # In the first step this will over-write the values in the matrices 
            self._put_to_matrices(xold, Itheta)
            
            # Store the flows for the beginning of the step
            self.FlowStorage.append(Flows_temporary)
            
            if Itheta >= 0.98*self.T.shape[1]:
                debug_plots(self)
                raise ValueError('98% of the maximum length of self.T reached, stopping calculation')
            
            t0+=h
            Itheta+=1
            xold = xnew
            
            #The error is already below the threshold
            if (max_error < eps_allowed and disableAdaptive == False):
#                print 'upsizing',h,
                #Take a bigger step next time, since eps_allowed>max_error
                h *= step_relax*(eps_allowed/max_error)**(0.2)
#                print h, eps_allowed, max_error
        
        #Store crank angle at the last index (first run at Itheta=0)
        self.t[Itheta] = t0
        # Re-evaluate derivs at the starting value for the step in order 
        # to use the correct values in the storage containers
        self.derivs(t0, xold)
        # Store the values for volumes and state vars in the matrices
        self._put_to_matrices(xold, Itheta)
        # Store the flows for the end
        self.FlowStorage.append(self.Flows.get_deepcopy())

        if sorted(self.stateVariables) == ['D','T']:
            self.CVs.updateStates('T',xnew[0:self.CVs.Nexist],'D',xnew[self.CVs.Nexist:2*self.CVs.Nexist])
        elif sorted(self.stateVariables) == ['M','T']:
            self.CVs.updateStates('T',xnew[0:self.CVs.Nexist],'D',xnew[self.CVs.Nexist:2*self.CVs.Nexist]/self.core.V)
        else:
            raise NotImplementedError
        
        # last index is Itheta, number of steps is Itheta+1
        print 'Itheta steps taken', Itheta+1
        self.Itheta=Itheta
        self.Ntheta = Itheta+1
        self.post_cycle()
        
    def calc_boundary_work(self):
        """
        This method calculates the boundary work rate using a trapezoidal 
        integration of
        
        .. math::
        
            \\dot W_{pv} = -\int p\\frac{dV}{d\\theta}\\frac{\\omega}{2\\pi} d\\theta
            
        for all the control volumes and sets the parameter ``self.Wdot_pv`` with 
        the result.
        
        The units of the boundary work are kW.
        """
        
        def Wdot_one_CV(CVindex):
            """ calculate the p-v work for one CV """
            
            x0_raw = self.V[CVindex, 0:self.Ntheta]
            y0_raw = self.p[CVindex, 0:self.Ntheta]
            
            x0, y0 = delimit_vector(x0_raw, y0_raw)
            
            summer=0.0
            for x_chunk, y_chunk in zip(x0, y0):
                summer += -trapz(y_chunk, x_chunk)*self.omega/(2*pi)
           
            return summer
            
        self.Wdot_pv = 0.0
        for CVindex in range(self.p.shape[0]):
            self.Wdot_pv+=Wdot_one_CV(CVindex)
        
    def post_cycle(self):
        """
        This stuff all happens at the end of the cycle.  It is a private method 
        not meant to be called externally
        
        The following things are done:
        #. The boundary work is calculated
        #. The flows are post-processed
        #. The heat transfer is post-processed
        #. The mass flow rate is calculated
        #. The volumetric efficiency is calculated
        #. The adiabatic efficiency is calculated
        #. The isentropic power is calculated
        #. The power input is calculated
        """
            
        self.calc_boundary_work()
            
        self._postprocess_flows()
        self._postprocess_HT()
        
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

        # Can't use intermediate temperature because the state might be two-phase
        # for some conditions and you are better off just calculating the enthalpy
        # directly
        h2s = Props('H','P',outletState.p,'S',s1,outletState.Fluid)
        
        self.eta_a = (h2s-h1)/(h2-h1)
        self.Wdot_i = self.mdot*(h2s-h1)
        # self.Qamb is positive if heat is being added to the lumped mass
        self.Wdot = self.mdot*(h2-h1)-self.Qamb
    
    def _check_cycle_abort(self, index, I = 100):
        """
        This function will check whether an abort has been requested every 
        ``I`` steps of the solver throughout the rotation
        
        Meant for calling by cycle_RK45, cycle_SimpleEuler, cycle_Heun, etc.
        
        Primarily this is useful for use with the GUI, where the GUI can pass
        an abort command to the model
        
        Parameters
        ----------
        index : int
            The index of the step
        I : int, optional
            Check abort at this interval
        
        """
        # % is the symbol for modulus in python
        if index % I == 0 and self.Abort():
            self._want_abort = True
            return True
        
    def check_abort(self):
        """
        A callback for use with the graphical user interface to force solver to quit
        
        It will check the Scroll.pipe_abort pipe for a ``True`` value, and if it
        finds one, it will set the Scroll._want_abort value to ``True`` which 
        will be read by the main execution thread
        
        Once self._want_abort is ``True``, it will stay latched True until the 
        run is terminated
        """
        #If you received an abort request, set a flag in the simulation
        if self.pipe_abort.poll() and self.pipe_abort.recv():
            print 'received an abort request'
            self._want_abort = True
            
        return self._want_abort
        
    def precond_solve(self,**kwargs):
        """
        A preconditioned solve
        
        What happens in this function is that one cycle is run using the very 
        rough first estimate of outlet state.  This cycle is then used to update
        the guess value for the outlet state.
        
        All keyword arguments are passed on to the PDSimCore.solve
        function
        """
        
        # We are going to over-write OldCycle in the preconditioner
        # so we cache its current value.  If not provided, use False
        OneCycle_oldval = kwargs.get('OneCycle', False)
        
        #Run solve with OneCycle turned on
        kwargs['OneCycle'] = True
        solve_output = self.solve(**kwargs)
        
        #If abort has already been called, stop and don't keep going
        if not self._want_abort:
            
            #Iterate over the tubes in order to find the inlet and outlet states
            for key, State in self.Tubes.Nodes.iteritems():
                if key == self.key_inlet:
                     IS = State
                if key == self.key_outlet:
                     OS = State
                     
            # Using Wdot_pv (boundary work), calculate a good guess for the 
            # outlet enthalpy by assuming the compressor is adiabatic
            h2 = IS.h + self.Wdot_pv/self.mdot
            
            # Find temperature as a function of enthalpy
            T2 = Props('T','H',h2,'P',OS.p,OS.Fluid)
            
            # Now run using the old value of OldCycle if OneCycle was 
            # not equal to True originally
            if not OneCycle_oldval:
                kwargs['OneCycle'] = OneCycle_oldval
                kwargs['x0']  = [T2, 0.5*T2+0.5*self.Tamb] #guesses for temperatures of discharge and lump
                kwargs['reset_initial_state'] = True
                self.solve(**kwargs)
       
    def connect_callbacks(self,
                          step_callback = None,
                          heat_transfer_callback = None,
                          lumps_energy_balance_callback = None,
                          endcycle_callback = None
                          ):
        """ 
        Connect up the callbacks for the simulation
        
        The callbacks must either be unbound methods or methods of a class derived from PDSimCore
        
        No keyword arguments are supported to be passed to the callbacks.  The 
        callback is probably a bound method of a PDSimCore instance, in which 
        case you have access to all the data in the class anyway
        
        Parameters
        ----------
        step_callback : function, or :class:`StepCallback <PDSim.core.callbacks.StepCallback>` subclass
            
            If a function is provided, it must have the call signature::
            
                disable_adaptive,h = step_callback(double t, double h, int i)
                
            where ``h`` is the step size that the adaptive solver wants to use, ``t`` is the current value of the independent variable, and ``i`` is the index in the container variables.  The return value ``disableAdaptive`` is a boolean value that describes whether the adaptive method should be turned off for this step ( ``False`` : use the adaptive method), and ``h`` is the step size you want to use.  If you don't want to disable the adaptive method and use the given step size, just::
                
                return False,h
            
            in your code.
        
        heat_transfer_callback : function, or :class:`HeatTransferCallback <PDSim.core.callbacks.HeatTransferCallback>` subclass
            
            If a function is provided, the heat_transfer_callback function must have the call signature::
            
                Q = heat_transfer_callback(double t)
            
            It should return an :class:`arraym <PDSim.misc.datatypes.arraym>` instance 
            with the same length as the number of CV in existence.  
            The entry in the :class:`arraym <PDSim.misc.datatypes.arraym>` is 
            positive if the heat transfer is TO the fluid in the CV in order 
            to maintain the sign convention that energy (or mass) input is 
            positive.
            
        lumps_energy_balance_callback : function, or :class:`LumpsEnergyBalanceCallback <PDSim.core.callbacks.LumpsEnergyBalanceCallback>` subclass
            
            If a function is provided, the lumps_energy_balance_callback 
            function must have the call signature::
            
                r = lumps_energy_balance_callback()
            
            It should return an :class:`arraym <PDSim.misc.datatypes.arraym>` 
            instance with the same length as the number of lumps.  The entry in 
            ``r`` is the value of the energy balance.  It will be driven to zero 
            by the solver
            
        """
        
        if step_callback is None:
            #No callback is provided, don't do anything
            pass
        elif isinstance(step_callback, PDSim.core.callbacks.StepCallback):
            #If the cythonized step callback is provided, hold onto it
            self.callbacks.step_callback = step_callback
        #Otherwise, wrap the desired callback if it has the right signature
        else:
            #Check the functional call
            callargs = inspect.getcallargs(step_callback, 0.0, 1e-10, 0)
            
            # Either a non-bound method is provided, or bound method is provided, in which case you get self,t,h,i as the values
            # t is a subclass of float, h is a subclass of float, is a subclass of int, and self is subclass of PDSimCore
            if not all([isinstance(arg,(float,int,PDSimCore)) for arg in callargs.values()]):
                sig_ok = False
            else:
                if len(callargs) in [3,4]:
                    sig_ok = True
                else:
                    sig_ok = False
            
            if step_callback is not None and sig_ok:
                self.callbacks.step_callback = PDSim.core.callbacks.WrappedStepCallback(self, step_callback)
            else:
                raise ValueError("step_callback is not possible to be wrapped - neither a subclass of StepCallback nor acceptable function signature")
        
        if heat_transfer_callback is None:
            #No callback is provided, don't do anything
            pass
        elif isinstance(heat_transfer_callback, PDSim.core.callbacks.HeatTransferCallback):
            #If the cythonized heat transfer callback is provided, hold a pointer to it
            self.callbacks.heat_transfer_callback = heat_transfer_callback
        else:
            callargs = inspect.getcallargs(heat_transfer_callback, 0.0)
            # Either a non-bound method is provided, or bound method is provided, in which case you get self,t as the values
            # t is a subclass of float, and self is subclass of PDSimCore
            if not all([isinstance(arg,(float,int,PDSimCore)) for arg in callargs.values()]):
                sig_ok = False
            else:
                if len(callargs) in [1,2]:
                    sig_ok = True
                else:
                    sig_ok = False
            
            #Otherwise, wrap the desired callback if it has the right signature
            if heat_transfer_callback is not None and sig_ok:
                self.callbacks.heat_transfer_callback = PDSim.core.callbacks.WrappedHeatTransferCallback(self, heat_transfer_callback)
            else:
                raise ValueError("heat_transfer_callback is not possible to be wrapped - neither a subclass of HeatTransferCallback nor an acceptable function")
        
        if lumps_energy_balance_callback is None:
            #No callback is provided, don't do anything
            pass
        elif isinstance(lumps_energy_balance_callback, PDSim.core.callbacks.LumpsEnergyBalanceCallback):
            #If the cythonized lump energy balance callback is provided, hold onto it
            self.callbacks.lumps_energy_balance_callback = lumps_energy_balance_callback
        #Otherwise, wrap the desired callback if it has the right signature
        else:
            callargs = inspect.getcallargs(lumps_energy_balance_callback)
            # Either a non-bound method is provided, or bound method is provided, in which case you get self,t as the values
            # t is a subclass of float, and self is subclass of PDSimCore
            sig_ok = len(callargs) == 0 or (len(callargs) == 1 and isinstance(callargs.values()[0],PDSimCore))
            
            if lumps_energy_balance_callback is not None and sig_ok:  #Do functional introspection here where the ``True`` is
                self.callbacks.lumps_energy_balance_callback = PDSim.core.callbacks.WrappedLumpsEnergyBalanceCallback(self, lumps_energy_balance_callback)
            else:
                raise ValueError("lump_energy_balance_callback is not possible to be wrapped - neither a subclass of LumpsEnergyBalanceCallback nor an acceptable function")
        
        self.callbacks.endcycle_callback = endcycle_callback
        
    def solve(self,
              key_inlet = None,
              key_outlet = None,
              solver_method = 'Euler',
              OneCycle = False,
              Abort = None,
              pipe_abort = None,
              UseNR = True,
              alpha = 0.5,
              plot_every_cycle = False,
              x0 = None,
              reset_initial_state = False,
              LS_start = 18,
              **kwargs):
        """
        This is the driving function for the PDSim model.  It can be extended through the 
        use of the callback functions
        
        It is highly recommended to call this function using keyword arguments like::
        
            solve(key_inlet = 'inlet.1', 
                  key_outlet = 'outlet.1', ....)
        
        Parameters
        ----------
        key_inlet : string
            The key for the flow node that represents the upstream quasi-steady point
        key_outlet : string
            The key for the flow node that represents the upstream quasi-steady point
        solver_method : string
        OneCycle : boolean
            If ``True``, stop after just one rotation.  Useful primarily for 
            debugging purposes
        Abort : function
            A function that may be called to determine whether to stop running.  
            If calling Abort() returns ``True``, stop running 
        pipe_abort : 
        UseNR : boolean
            If ``True``, use a multi-dimensional solver to determine the initial state of the state variables for each control volume
        alpha : float
            Use a range of ``(1-alpha)*dx, (1+alpha)*dx`` for line search if needed
        plot_every_cycle : boolean
            If ``True``, make the plots after every cycle (primarily for debug purposes)
        x0 : arraym
            The starting values for the solver that modifies the discharge temperature and lump temperatures
        reset_initial_state : boolean
            If ``True``, use the stored initial state from the previous call to ``solve`` as the starting value for the thermodynamic values for the control volumes
        LS_start : int
            Number of conventional steps to be taken when not using newton-raphson prior to entering into a line search
            
        step_callback : function
            DEPRECATED! Should be passed to the connect_callbacks() function before running precond_solve() or solve()
        endcycle_callback : function
            DEPRECATED! Should be passed to the connect_callbacks() function before running precond_solve() or solve()
        heat_transfer_callback : function
            DEPRECATED! Should be passed to the connect_callbacks() function before running precond_solve() or solve()
        lump_energy_balance_callback : function
            DEPRECATED! Should be passed to the connect_callbacks() function before running precond_solve() or solve()
        valves_callback : function
            DEPRECATED! Should be passed to the connect_callbacks() function before running precond_solve() or solve()
        """
        if any(cb in kwargs for cb in ['step_callback','endcycle_callback','heat_transfer_callback','lump_energy_balance_callback','valves_callback']):
            raise NotImplementedError('callback functions are no longer passed to solve() function, rather they are passed to connect_callbacks() function prior to calling solve()')
        
        #Connect functions that have been serialized by saving the function name as a string
        self.connect_flow_functions()
        
        #Both inlet and outlet keys must be connected to invariant nodes - 
        # that is they must be part of the tubes which are all quasi-steady
        if not key_inlet == None and not key_inlet in self.Tubes.Nodes:
            raise KeyError('key_inlet must be a Tube node')
        if not key_outlet == None and not key_outlet in self.Tubes.Nodes:
            raise KeyError('key_outlet must be a Tube node')
        
        self.key_inlet = key_inlet
        self.key_outlet = key_outlet
        
        from time import clock
        t1=clock()
        
        # Set up a pipe for accepting a value of True which will abort the run
        # Used from the GUI to kill process from the top-level thread
        self.pipe_abort = pipe_abort
        
        # If a function called pre_solve is provided, call it with no input arguments
        if hasattr(self,'pre_solve'):
            self.pre_solve()
            
        # This runs before the model starts at all
        self.pre_run()
        
        # Check which method is used to do aborting
        if Abort is None and pipe_abort is not None:
            # Use the pipe_abort pipe to look at the abort pipe to see whether 
            # to quit 
            self.Abort = self.check_abort
        elif Abort is None and pipe_abort is None:
            #Disable the ability to abort, always don't abort
            self.Abort = lambda : False
        elif Abort is not None and pipe_abort is None:
            self.Abort = Abort
        else:
            raise ValueError('Only one of Abort and pipe_abort may be provided')
        
        # If you want to reset the initial state, use the values that were
        # cached in the xstate_init array
        if reset_initial_state is not None and reset_initial_state:
            self.reset_initial_state()
            self.pre_cycle(self.xstate_init)
        else:
            # (2) Run a cycle with the given values for the temperatures
            self.pre_cycle()
            #x_state only includes the values for the chambers, the valves start closed
            #Since indexed, yields a copy
            self.x_state = self._get_from_matrices(0)[0:len(self.stateVariables)*self.CVs.Nexist]
        
        if x0 is None:
            x0 = [self.Tubes.Nodes[key_outlet].T, self.Tubes.Nodes[key_outlet].T]
                
        def OBJECTIVE_ENERGY_BALANCE(Td_Tlumps):
            """
            Td_Tlumps: arraym instance
                Contains the discharge temperature followed by the temperatures of each lumped mass
            """
            print Td_Tlumps,'Td,Tlumps inputs'
            # Td_Tlumps is a list (or or np.ndarray)
            Td_Tlumps = list(Td_Tlumps)
            # Consume the first element as the discharge temp 
            self.Td = float(Td_Tlumps.pop(0))
            # The rest are the lumps in order
            self.Tlumps = Td_Tlumps
            
            # (0). Update the discharge temperature
            p = self.Tubes.Nodes[key_outlet].p
            self.Tubes.Nodes[key_outlet].update({'T':self.Td,'P':p})
            
            def OBJECTIVE_CYCLE(x_state):
                """
                x_state: arraym instance
                    Contains the state variables for all the control volumes in existence, as well as the valve values
                """
                x_state = arraym(x_state)
                
                # (1). First, run all the tubes
                for tube in self.Tubes:
                    tube.TubeFcn(tube)
                    
                # Call update_existence to save the enthalpies for the tubes 
                self.update_existence()
                
                #The first time this function is run, save the initial state
                # and the existence of the CV, as well as the valve positions
                if self.xstate_init is None:
                    self.xstate_init = x_state
                    self.exists_CV_init = self.CVs.exists_keys
                try:
                    t1 = clock()
                    if solver_method == 'Euler':
                        #Default to 7000 steps if not provided
                        N = getattr(self,'EulerN', 7000)
                        aborted = self.cycle_SimpleEuler(N,x_state)
                    elif solver_method == 'Heun':
                        #Default to 7000 steps if not provided
                        N = getattr(self,'HeunN', 7000)
                        aborted = self.cycle_Heun(N,x_state)
                    elif solver_method == 'RK45':
                        #Default to tolerance of 1e-8 if not provided
                        eps_allowed = getattr(self,'RK45_eps', 1e-8)
                        aborted = self.cycle_RK45(x_state,eps_allowed = eps_allowed)
                    else:
                        raise AttributeError('solver_method should be one of RK45, Euler, or Heun')
                except ValueError:
                    debug_plots(self)
                    raise
                
                t2 = clock()
                print 'Elapsed time for cycle is {0:g} s'.format(t2-t1)
                
                #Quit if you have aborted in one of the cycle solvers
                if aborted == 'abort':
                    return None
                    
                mdot_out = self.FlowsProcessed.mean_mdot[key_outlet]
                mdot_in = self.FlowsProcessed.mean_mdot[key_inlet]
                print 'Mass flow difference',(mdot_out+mdot_in)/mdot_out*100,' %'
                
                # We need to find the key at the inlet to the outlet tube.
                for Tube in self.Tubes:
                    if Tube.key1 == key_outlet:
                        key_outtube_inlet = Tube.key2
                        break
                    elif Tube.key2 == key_outlet:
                        key_outtube_inlet = Tube.key1
                        break
                    
                # This is the so-called hd' state at the outlet of the pump set
                h_outlet = (self.FlowsProcessed.integrated_mdoth[key_outtube_inlet]
                            /self.FlowsProcessed.integrated_mdot[key_outtube_inlet])
                # It should be equal to the enthalpy of the fluid at the inlet
                # to the outlet tube at the current Td value
                h_outlet_Tube = self.Tubes.Nodes[key_outtube_inlet].h
                # Residual is the difference of these two terms
                # We put it in kW by multiplying by flow rate
                self.resid_Td = mdot_out * (h_outlet_Tube - h_outlet)
                
                if plot_every_cycle:
                    debug_plots(self)
                
                #If the abort function returns true, quit this loop
                if self.Abort() or OneCycle:
                    print 'Quitting OBJECTIVE_CYCLE loop in core.solve'
                    return None #Stop
                        
                if self.callbacks.endcycle_callback is None:
                    return None #Stop
                else:
                    #endcycle_callback returns the errors and new initial state for the solver
                    errors, x_state_new = self.callbacks.endcycle_callback()
                    self.x_state = x_state_new.copy() #Make a copy
                    return errors
                
                ##################################
                ## End OBJECTIVE_CYCLE function ##
                ##################################
            
            if UseNR:
                self.x_state = Broyden(OBJECTIVE_CYCLE,
                                       self.x_state, 
                                       Nfd = 1, 
                                       dx = 0.01*arraym(self.x_state), 
                                       itermax = 50, 
                                       ytol = 1e-4)
                if self.x_state[0] is None:
                    return None
            else:
                diff_old = None
                x_state_prior = None
                
                init_state_counter = 0
                while True:
                    # This block runs the first time through 
                    # (when old state is not defined)
                    # This means taking a full step
                    if x_state_prior is None or init_state_counter < LS_start:
                        x_state_prior = self.x_state.copy()
                        errors = OBJECTIVE_CYCLE(x_state_prior)
                        if errors is None:
                            break
                        else:
                            errors *= 100
                            x_state_new = self.x_state.copy()
                            RSSE_prior = np.sqrt(np.sum(np.power(errors, 2)))
                    try:  
                        
                        if init_state_counter >= LS_start:
                            dx = x_state_new - x_state_prior
                            print textwrap.dedent(
                            """
                            ****************************
                            Starting a rough line search
                            ****************************
                            """
                            )
                            # If you are not getting good convergence, 
                            # use a rough linear search
                            _w, _error = [0.0], [RSSE_prior]
                            for w in [0.5, 1.0, 1.5, 2.0, 2.5]:
                                x = x_state_prior + w*dx
                                errors  = OBJECTIVE_CYCLE(x)*100
                                RSSE = np.sqrt(np.sum(np.power(errors, 2)))
                                _w.append(w)
                                _error.append(RSSE)
                                print 'w',w,'RSSE',RSSE
#                            f = plt.figure()
#                            ax=f.add_axes((0.15,0.15,0.8,0.8))
#                            ax.plot(_w, _error)
#                            plt.show()
                            i = _error.index(min(_error))
                            x_state_new = x_state_prior + _w[i]*dx
                            init_state_counter = 0
                            
                        else:
                            if errors is not None:
                                if not len(x_state_new) == len(x_state_prior):
                                    raise IndexError
                                diff = (x_state_prior - x_state_new) / x_state_prior * 100
                                diff_abs = [abs(dx) for dx in diff]
                                RSSE = np.sqrt(np.sum(np.power(errors, 2)))
                                print 'Cycle #',init_state_counter,'RSSE',RSSE

                    except IndexError:
                        # You will get an IndexError if the length of the x_state
                        # list changes due to different CV being in existence at
                        # the beginning and end of the rotation.  Can be caused
                        # by a volume ratio that puts the discharge angle right
                        # before the end of the rotation
                        continue

                    # If the attribute eps_cycle is set, use that to determine 
                    # the convergence criterion, otherwise, use the default
                    if RSSE < getattr(self,'eps_cycle',1e-3):
                        break
                    else:
                        self.x_state = x_state_new
                    
                    #Increment the counter for the inner loop
                    init_state_counter += 1
            
            #If the abort function returns true, quit this loop
            if self.Abort() or OneCycle:
                print 'Quitting OBJECTIVE function in core.solve'
#                if self.callbacks.lumps_energy_balance_callback is not None:
#                    self.callbacks.lumps_energy_balance_callback()
                return None
                    
            # (3) After convergence of the inner loop, check the energy balance on the lumps
            if self.callbacks.lumps_energy_balance_callback is not None:
                resid_HT = self.callbacks.lumps_energy_balance_callback()
                if not isinstance(resid_HT,arraym):
                    resid_HT = arraym(resid_HT)
                
            resids = [self.resid_Td]
            resids.extend(resid_HT) #Need to extend because resid_HT is an arraym
            print resids,'resids'
            return resids
        
            ##################################
            ## End OBJECTIVE_ENERGY_BALANCE ##
            ##################################
        
        x_soln = Broyden(OBJECTIVE_ENERGY_BALANCE,x0, dx=1.0, ytol=0.01, itermax=30)
        print 'Solution is', x_soln,'Td, Tlumps'
        
        if not self.Abort(): 
            self.post_solve()
            
        if hasattr(self,'resid_Td') and hasattr(self,'x_state'):
            del self.resid_Td, self.x_state
        
        # Save copies of the inlet and outlet states at the root of the HDF5 file
        self.inlet_state = self.Tubes.Nodes[key_inlet]    
        self.outlet_state = self.Tubes.Nodes[key_outlet]
            
        #Save the elapsed time for simulation
        self.elapsed_time = clock()-t1
        
    def get_prune_keys(self):
        """
        Remove some elements when the simulation finishes that are not 
        very useful and/or are very large when stored to file
        
        Returns
        -------
        prune_key_list: list
            A list of HDF5 keys that are to be removed from the HDF5 file.
        """
        
        return ['/CVs/CVs',
                '/CVs/Nodes',
                '/CVs/T',
                '/CVs/cp',
                '/CVs/cv',
                '/CVs/dpdT',
                '/CVs/exists_CV',
                '/CVs/exists_indices',
                '/CVs/exists_keys',
                '/CVs/h',
                '/CVs/p',
                '/CVs/rho',
                '/callbacks',
                '/core',
                '/steps',
                '/theta',
                '/Flows'
                ]
        
    def attach_HDF5_annotations(self, fName):
        """
        In this function, annotations can be attached to each HDF5 field
        
        Parameters
        ----------
        fName : string
            The file name for the HDF5 file that is to be used
        """ 
        attrs_dict = {
                '/t':'The array of the independent variable in the solution, either time or crank angle [rad or s]',
                '/m':'The NCV x Nt matrix with the mass in each control volume [kg]',
                '/T':'The NCV x Nt matrix with the temperature in each control volume [K]',
                '/V':'The NCV x Nt matrix with the volume in each control volume [m^3]',
                '/dV':'The NCV x Nt matrix with the derivative of volume w.r.t. crank angle in each control volume [m^3/radian]',
                '/h':'The NCV x Nt matrix with the enthalpy in each control volume [kJ/kg]',
                '/p':'The NCV x Nt matrix with the pressure in each control volume [kPa]',
                '/rho':'The NCV x Nt matrix with the density in each control volume [kg/m^3]',
                '/Q':'The NCV x Nt matrix with the heat transfer TO the gas in each control volume [kW]',
                '/xL':'The NCV x Nt matrix with the oil mass fraction in each control volume [-]',
                '/A_shell':'The shell area of the machine [m^2]',
                '/h_shell':'The heat transfer coefficient between the shell and the ambient [kW/m^2/K]',
                '/key_inlet':'The key for the inlet node',
                '/key_outlet':'The key for the outlet node',
                '/elapsed_time':'The elapsed time for the simulation run [s]',
                '/eta_a': 'Adiabatic efficiency [-]',
                '/eta_oi':'Overall isentropic efficiency [-]',
                '/eta_v':'Volumetric efficiency [-]',
                '/Qamb':'Rate of heat transfer from the machine TO the ambient [kW]',
                '/RK45_eps':'Step error tolerance for Runge-Kutta method [varied]',
                '/Tamb':'Ambient temperature [K]',
                '/Wdot_pv':'Mechanical power calculated as the integral of -pdV [kW]',
                '/Wdot_electrical':'Electrical power of the machine [kW]',
                '/Wdot_forces':'Mechanical power calculated from the mechanical analysis [kW]',
                '/motor/eta_motor':'Motor efficiency [-]',
                '/motor/losses':'Losses generated in the motor [kW]',
                '/motor/suction_fraction':'Fraction of the motor losses that are added to the suction gas [-]',
                '/motor/type':'The model used to simulate the motor',
                '/run_index':'A unique identifier for runs in a batch'
                }
        
        hf = h5py.File(fName,'a')
        
        for k, v in attrs_dict.iteritems():
            dataset = hf.get(k)
            if dataset is None:
                print 'bad key',k
            else:
                dataset.attrs['note'] = v
        hf.close()
        
    def post_solve(self):
        """
        Do some post-processing to calculate flow rates, efficiencies, etc.  
        """      
        
        #Resize all the matrices to keep only the real data
        print 'Ntheta is', self.Ntheta
        self.t = self.t[  0:self.Ntheta]
        self.T = self.T[:,0:self.Ntheta]
        self.p = self.p[:,0:self.Ntheta]
        self.Q = self.Q[:,0:self.Ntheta]
        self.m = self.m[:,0:self.Ntheta]
        self.rho = self.rho[:,0:self.Ntheta]
        self.V = self.V[:,0:self.Ntheta]
        self.dV = self.dV[:,0:self.Ntheta]
        self.h = self.p[:,0:self.Ntheta]
        self.xL = self.p[:,0:self.Ntheta]
        self.xValves = self.xValves[:,0:self.Ntheta]
        
        
        
        print 'mdot*(h2-h1),P-v,Qamb', self.Wdot, self.Wdot_pv, self.Qamb
        print 'Mass flow rate is',self.mdot*1000,'g/s'
        print 'Volumetric efficiency is',self.eta_v*100,'%'
        print 'Adiabatic efficiency is',self.eta_a*100,'%'
        
    def derivs(self,theta,x):
        """
        Evaluate the derivatives of the state variables
        
        derivs() is an internal function that should (probably) not actually be called by
        any user-provided code, but is documented here for completeness.
        
        Parameters
        ----------
        theta : float
            The value of the independent variable
        x : ``arraym`` instance
            The array of the state variables (plus valve parameters)
        
        Returns
        -------
        dfdt : ``arraym`` instance
        
        """

        # Updates the state, calculates the volumes, prepares all the things needed for derivatives
        self.core.properties_and_volumes(self.CVs.exists_CV, theta, STATE_VARS_TM, x)
        
        #Join the enthalpies of the CV in existence and the tubes
        harray = self.core.h.copy()
        harray.extend(self.Tubes.get_h())
        
        # Calculate the flows and sum up all the terms
        self.core.calculate_flows(self.Flows, harray)
        
        # Calculate the heat transfer terms if provided
        if self.callbacks.heat_transfer_callback is not None:
            self.core.Q = arraym(self.callbacks.heat_transfer_callback(theta))
            if not len(self.core.Q) == self.CVs.Nexist:
                raise ValueError('Length of Q is not equal to length of number of CV in existence')
        else:
            self.core.Q = empty_arraym(self.CVs.Nexist)
            
        # Calculate the derivative terms and set the derivative of the state vector
        self.core.calculate_derivs(self.omega, False)
        
        #Liquid not yet supported
        if self.__hasLiquid__:
            raise NotImplementedError()
        
        #Add the derivatives for the valves
        if self.__hasValves__:
            # 
            offset = len(self.stateVariables)*self.CVs.Nexist
            for i,valve in enumerate(self.Valves):
                #Get the values from the input array for this valve
                xvalve = x[offset+i*2:offset+2+i*2]
                #Set the values in the valve class
                valve.set_xv(xvalve)
                # Get the derivatives of position and derivative of velocity
                self.core.property_derivs.extend(valve.derivs(self))
        
        return self.core.property_derivs
        
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
        
    def IsentropicNozzleFM(self,FlowPath,A,**kwargs):
        """
        A generic isentropic nozzle flow model wrapper
        
        Parameters
        ----------
        FlowPath : FlowPath instance
            A fully-instantiated flow path model
        A : float
            throat area for isentropic nozzle model [:math:`m^2`]
            
        Returns
        -------
        mdot : float
            The mass flow through the flow path [kg/s]
        """

        try:
            mdot = flow_models.IsentropicNozzle(A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
        
    def IsentropicNozzleFMSafe(self,FlowPath,A,DP_floor,**kwargs):
        """
        A generic isentropic nozzle flow model wrapper with the added consideration
        that if the pressure drop is below the floor value, there is no flow.
        This was added to handle the case of the injection line where there is
        no flow out of the injection which greatly increases the numerical 
        stiffness 
        
        Parameters
        ----------
        FlowPath : FlowPath instance
            A fully-instantiated flow path model
        A : float
            throat area for isentropic nozzle model [:math:`m^2`]
        DP_floor: float
            The minimum pressure drop [kPa]
            
        Returns
        -------
        mdot : float
            The mass flow through the flow path [kg/s]
        """

        try:
            if FlowPath.State_up.p-FlowPath.State_down.p > DP_floor:
                mdot = flow_models.IsentropicNozzle(A,
                                                    FlowPath.State_up,
                                                    FlowPath.State_down)
                return mdot
            else:
                
                return 0.0
        except ZeroDivisionError:
            return 0.0

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
        assert self.Ntheta - 1 == self.Itheta
        #old and new CV keys
        LHS,RHS=[],[]
        errorT,error_rho,error_mass,newT,new_rho,new_mass,oldT,old_rho,old_mass={},{},{},{},{},{},{},{},{}
        
        for key in self.CVs.exists_keys:
            # Get the 'becomes' field.  If a list, parse each fork of list. If a single key convert 
            # into a list so you can use the same code below 
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
                #What the state variables were at the start of the rotation
                oldT[newkey]=self.T[Inew, 0]
                old_rho[newkey]=self.rho[Inew, 0]
                #What they are at the end of the rotation
                newT[newkey]=self.T[Iold,self.Itheta]
                new_rho[newkey]=self.rho[Iold,self.Itheta]
                
                errorT[newkey]=(oldT[newkey]-newT[newkey])/newT[newkey]
                error_rho[newkey]=(old_rho[newkey]-new_rho[newkey])/new_rho[newkey]
                #Update the list of keys for setting the exist flags
                LHS.append(key)
                RHS.append(newkey)
        
        error_T_list = [errorT[key] for key in self.CVs.keys if key in newT]
        error_rho_list = [error_rho[key] for key in self.CVs.keys if key in new_rho]
        
        new_T_list = [newT[key] for key in self.CVs.keys if key in newT]
        new_rho_list = [new_rho[key] for key in self.CVs.keys if key in new_rho]
        
        #Reset the exist flags for the CV - this should handle all the possibilities
        #Turn off the LHS CV
        for key in LHS:
            self.CVs[key].exists=False
        #Turn on the RHS CV
        for key in RHS:
            self.CVs[key].exists=True
            
        self.update_existence()

        # Error values are based on density and temperature independent of 
        # selection of state variables 
        error_list = []
        for var in ['T','D']:
            if var == 'T':
                error_list += error_T_list
            elif var == 'D':
                error_list += error_rho_list
            elif var == 'M':
                error_list += error_mass_list
            else:
                raise KeyError
            
        # Calculate the volumes at the beginning of the next rotation
        self.core.just_volumes(self.CVs.exists_CV, 0)
        V = {key:V for key,V in zip(self.CVs.exists_keys,self.core.V)}
        new_mass_list = [new_rho[key]*V[key] for key in self.CVs.exists_keys]
        
        new_list = []
        for var in self.stateVariables:
            if var == 'T':
                new_list += new_T_list
            elif var == 'D':
                new_list += new_rho_list
            elif var == 'M':
                new_list += new_mass_list
            else:
                raise KeyError
    
        return arraym(error_list), arraym(new_list)
    
    def connect_flow_functions(self):
        """
        Reconnect function pointers
        
        For pickling purposes, it can sometimes be useful to just give the name
        of the function relative to the PDSimCore (or derived class).  If the 
        function is just a string, reconnect it to the function in the PDSimCore 
        instance 
        """
        for Flow in self.Flows:
            if hasattr(Flow.MdotFcn, 'Function'):
                if isinstance(Flow.MdotFcn.Function, basestring):
                    if hasattr(self,Flow.MdotFcn.Function):
                        Flow.MdotFcn.Function = getattr(self, Flow.MdotFcn.Function)
                    else:   
                        raise AttributeError('The name of the function ['+Flow.MdotFcn.Function+']is not found in the PDSimCore derived class instance')
            
if __name__=='__main__':
    PC = PDSimCore()
    PC.attach_HDF5_annotations('runa.h5')
    print 'This is the base class that is inherited by other compressor types.  Running this file doesn\'t do anything'
