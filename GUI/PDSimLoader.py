# -*- coding: latin-1 -*-

from math import pi
from CoolProp import State as CPState
from PDSim.core.containers import ControlVolume,ControlVolumeCollection
from PDSim.flow.flow import FlowPath
from PDSim.flow.flow_models import ValveModel
from PDSim.core.core import Tube
from time import clock
import threading

def RecipBuilder(recip):
    
    recip.pre_solve()
    
    outletState=CPState.State(recip.inletState.Fluid,{'T':400,'P':recip.discharge_pressure})
    mdot_guess = recip.inletState.rho*recip.Vdisp()*recip.omega/(2*pi)
    
    #First add the control volumes.
    recip.add_CV( ControlVolume(key='A',
                               initialState=outletState.copy(),
                               VdVFcn=recip.V_dV,
                               becomes='A') )
    recip.add_CV( ControlVolume(key='shell',
                               initialState=recip.inletState.copy(),
                               VdVFcn=recip.V_shell,
                               becomes='shell') )
    
    recip.add_tube( Tube(key1='inlet.1',key2='inlet.2',L=0.03,ID=0.02,
                             mdot=mdot_guess, State1=recip.inletState.copy(),
                             fixed=1,TubeFcn=recip.TubeCode) )
    recip.add_tube( Tube(key1='outlet.1',key2='outlet.2',L=0.03,ID=0.02,
                             mdot=mdot_guess, State2=outletState.copy(),
                             fixed=2,TubeFcn=recip.TubeCode) )
    
    recip.add_flow(FlowPath(key1='shell', key2='inlet.2', MdotFcn=recip.Inlet))
    recip.add_flow(FlowPath(key1='shell', key2='A', MdotFcn=recip.Suction))
    recip.add_flow(FlowPath(key1='outlet.1', key2='A', MdotFcn=recip.Discharge))
    recip.add_flow(FlowPath(key1='shell', key2='A', MdotFcn=recip.PistonLeakage))
    
    #The suction valve parameters
    recip.suction_valve=ValveModel(
          d_valve=recip.valve_d,
          d_port=recip.d_suction,
          C_D=recip.valve_C_D,
          h_valve=recip.valve_h,
          a_valve=recip.valve_a,
          l_valve=recip.valve_l,
          rho_valve=recip.valve_rho,
          E=recip.valve_E,
          x_stopper=recip.valve_x_stopper,
          key_up='inlet.2',
          key_down='A'
          )
    recip.add_valve(recip.suction_valve)
    #The discharge valve parameters
    recip.discharge_valve=ValveModel(
          d_valve=recip.valve_d,
          d_port=recip.d_discharge,
          C_D=recip.valve_C_D,
          h_valve=recip.valve_h,
          a_valve=recip.valve_a,
          l_valve=recip.valve_l,
          rho_valve=recip.valve_rho,
          E=recip.valve_E,
          x_stopper=recip.valve_x_stopper,
          key_up='A',
          key_down='outlet.1'
          )
    recip.add_valve(recip.discharge_valve)
    
def ScrollBuilder(scroll):
    
    if hasattr(scroll,'pre_solve'):
        scroll.pre_solve()
    
    outletState=CPState.State(scroll.inletState.Fluid,{'T':400,'P':scroll.discharge_pressure})
    
    if callable(scroll.Vdisp):
        Vdisp = scroll.Vdisp()
    else:
        Vdisp = scroll.Vdisp
        
    mdot_guess = scroll.inletState.rho*Vdisp*scroll.omega/(2*pi)
    
    scroll.auto_add_CVs(scroll.inletState, outletState)
    scroll.auto_add_leakage(flankFunc = scroll.FlankLeakage, 
                            radialFunc = scroll.RadialLeakage)
    
    scroll.add_tube( Tube(key1='inlet.1',key2='inlet.2',L=0.03,ID=0.02,
                             mdot=mdot_guess, 
                             State1=scroll.inletState.copy(),
                             fixed=1,TubeFcn=scroll.TubeCode) )
    
    scroll.add_tube( Tube(key1='outlet.1',key2='outlet.2',L=0.03,ID=0.02,
                             mdot=mdot_guess, 
                             State2=outletState.copy(),
                             fixed=2,TubeFcn=scroll.TubeCode) )
    
    scroll.add_flow(FlowPath(key1='inlet.2',key2='sa',
                             MdotFcn=scroll.Inlet_sa))
    scroll.add_flow(FlowPath(key1='sa',key2='s1',MdotFcn=scroll.SA_S))
    scroll.add_flow(FlowPath(key1='sa',key2='s2',MdotFcn=scroll.SA_S))
    
    scroll.add_flow(FlowPath(key1='outlet.1',key2='ddd',MdotFcn=scroll.Discharge))
    scroll.add_flow(FlowPath(key1='outlet.1',key2='dd',MdotFcn=scroll.Discharge))
    scroll.add_flow(FlowPath(key1='d1',key2='dd',MdotFcn=scroll.D_to_DD))
    scroll.add_flow(FlowPath(key1='d2',key2='dd',MdotFcn=scroll.D_to_DD))
    
