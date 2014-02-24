""""
Example simulation for a reciprocating piston expander

(c) Ian Bell, 2013
"""

##########################################
###          Part 1. Imports           ###
##########################################

#Here we import some python packages
from __future__ import division
from math import pi, cos, sin
from time import clock
import os, sys

# If the following line is uncommented, python will try to use a local version
# of PDSim.  This is handy for debugging purposes.  Generally you want this line 
# commented out
#sys.path.insert(0, os.path.abspath('..'))

#Here we import the things from PDSim we need
from PDSim.flow.flow import FlowPath
from PDSim.flow import flow_models
from PDSim.misc.datatypes import empty_arraym
from PDSim.core.containers import ControlVolume, Tube
from PDSim.core.core import PDSimCore
from PDSim.plot.plots import debug_plots

#Imports from CoolProp (fluid property database)
from CoolProp import State
from CoolProp import CoolProp as CP
    
#################################################
###   Part 2. Declaration of PistonExpander   ###
#################################################

class PistonExpander(PDSimCore):
    
    #: Displacement of the cylinder above the dead volume [m^3]
    Vdisp = 20e-6
    
    #: Dead volume of the cylinder at TDC [m^3]
    Vdead = 3e-6
    
    #: Rotational speed [rad/s]
    omega = 377
    
    def __init__(self):
        #Initialize the base class that PistonExpander is derived from
        PDSimCore.__init__(self)
        
    def V_dV(self, theta):
        """
        The simplest volume relationship possible, given by dead volume 
        and displacement directly
        """
        V = self.Vdead+self.Vdisp/2*(1-cos(theta))
        dVdtheta = self.Vdisp/2*sin(theta)
        return V, dVdtheta
    
    def Suction(self, FlowPath):
        if 0 <= self.theta <= pi/4:
            FlowPath.A = pi*0.006**2/4*(1-cos(8*self.theta))/2
            mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                                FlowPath.State_up,
                                                FlowPath.State_down)
        else:
            FlowPath.A = 0.0
            mdot = 0
        return mdot
        
    def Discharge(self, FlowPath):
        if pi <= self.theta <= 7*pi/4:
            FlowPath.A = pi*0.006**2/4*(1-cos(4*self.theta))/2
            mdot = flow_models.IsentropicNozzle(FlowPath.A,
                                            FlowPath.State_up,
                                            FlowPath.State_down)
        else:
            FlowPath.A = 0.0
            mdot = 0
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
        
        More code (a proper heat transfer model) could be included here, like 
        in PDSim.recip.core
        """
        return empty_arraym(self.CVs.N)
        
    def mechanical_losses(self):
        """
        The mechanical losses in kW
        """
        return 0#self.Wdot_parasitic
    
    def ambient_heat_transfer(self):
        """
        The ambient heat transfer for the compressor in kW
        
        Returns a positive value if heat is added to the compressor from the 
        ambient
        """
        return self.h_shell*self.A_shell*(self.Tamb-self.Tlumps[0]) #[kW]
        
    def lump_energy_balance_callback(self):
        """
        A callback used in PDSimCore.solve to do the energy balance on the lump
        
        Note: we neglect heat transfer to the gas in the working chamber
        """
        #Mechanical losses are added to the lump
        self.Wdot_mechanical = self.mechanical_losses() #[kW]
        #Heat transfer between the shell and the ambient
        self.Qamb = self.ambient_heat_transfer() #[kW]
        return self.Wdot_mechanical + self.Qamb
    
    def step_callback(self,theta,h,Itheta):
        self.theta = theta
        return False, h
 
#######################################
###    Part 3. Execution of code    ###
#######################################
        
def Expander():

    expander = PistonExpander() #Instantiate the class
    
    Ref = 'Nitrogen'
    inletState = State.State(Ref,dict(T = 298.15, P = 501.325))
    outletState = State.State(Ref,dict(T = 200, P = inletState.p/10))
    mdot_guess = inletState.rho*expander.Vdisp*expander.omega/(2*pi)
    
    #First add the control volume.
    expander.add_CV(ControlVolume(key='A',
                                  initialState=inletState.copy(),
                                  VdVFcn=expander.V_dV,)
                    )
        
    #These are parameters needed for the ambient heat transfer model
    expander.h_shell = 0.010               #[kW/m2/K]
    expander.A_shell = pi*10*2*(0.0254**2) #[m2]
    expander.Tamb = 298                    #[K] 
    
    #Parameters for the mechanical losses model (simplified)
    expander.Wdot_parasitic = 0.01         #Parasitic losses [kW]
    
    #Add the inlet tube
    expander.add_tube(Tube(key1 = 'inlet.1',
                           key2 = 'inlet.2',
                           L = 0.03,
                           ID = 0.01,
                           mdot = mdot_guess, 
                           State1 = inletState.copy(),
                           fixed = 1,
                           TubeFcn = expander.TubeCode) 
                      )
    
    #Add the outlet tube
    expander.add_tube(Tube(key1 = 'outlet.1',
                           key2 = 'outlet.2',
                           L = 0.03,
                           ID = 0.01,
                           mdot = mdot_guess, 
                           State2 = outletState.copy(),
                           fixed = 2,
                           TubeFcn = expander.TubeCode) 
                      )
    
    #Add the flow paths that link flow nodes together
    expander.add_flow(FlowPath(key1='inlet.2',key2='A',MdotFcn=expander.Suction))
    expander.add_flow(FlowPath(key1='outlet.1',key2='A',MdotFcn=expander.Discharge))
    
    t1=clock()
    expander.EulerN = 4000
    expander.RK45_eps = 1e-10
    expander.connect_callbacks(step_callback = expander.step_callback,
                               endcycle_callback=expander.endcycle_callback, # Provided by PDSimCore
                               heat_transfer_callback=expander.heat_transfer_callback,
                               lumps_energy_balance_callback = expander.lump_energy_balance_callback)
    expander.solve(key_inlet='inlet.1',
                   key_outlet='outlet.2',
                   solver_method = 'Euler',
                   OneCycle = False,
                   UseNR = True,
                   plot_every_cycle = False
                   )
    print 'time taken',clock()-t1,'s'
    
    #debug_plots(expander)
    
if __name__=='__main__':    
    #If this file is run directly, this code will be run
    Expander()

    