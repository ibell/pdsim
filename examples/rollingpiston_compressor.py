"""
Domestic Hermetic Rolling Piston Compressor

References:
K.T. Ooi, "Design optimization of a rolling piston compressor for refrigerators", Applied Thermal Engineering, 25(2005), 813-829

K.T. Ooi, "Compressor Performance Comparison When Using R134a nad R1234yf as Working Fluids", International Compressor Engineering Conference at Purdue, 2012. paper 2146
"""
## This file is meant for experimentation with PDSim features and may not
## work for you.  Caveat emptor!!
#
from __future__ import division, print_function

import sys, os
from math import pi, cos, sin, asin,tan,sqrt

#PDSim imports
from PDSim.flow.flow import FlowPath
from PDSim.flow import flow_models
from PDSim.misc.datatypes import arraym,listm
from PDSim.core.containers import ControlVolume, Tube
from PDSim.core.core import PDSimCore,struct
from PDSim.plot.plots import debug_plots
#from PDSim.flow.flow_models import ValveModel
from PDSim.flow.ValveModels import ValveModel
from PDSim.core.motor import Motor
from PDSim.rolling.core_clean import RollingPiston

#CoolProp imports
from CoolProp import State,AbstractState
from CoolProp import CoolProp as CP

import numpy as np
import matplotlib.pyplot as plt
import timeit

def RollingPistonCompressor(Te = -10.6, Tc = 54.4, Ref = 'HEOS::R1234yf',HDF5file='rollingpiston_compressor.h5'):
    
    rolling=RollingPiston()                         #Instantiate the model
    #This runs if the module code is run directly
    
    rolling.e = 4.8e-3                              #Eccentricity, m
    rolling.Rc = 29.0e-3                            #Radius of cylinder, m
    rolling.Rc_outer = 35.0e-3                      #Radius of outer cylinder, m
    rolling.Rr = 23.4e-3                            #Radius roller, m
    rolling.Rr_i = 20e-3                            #Inner radius roller, m
    rolling.Re = 18.0e-3                            #Eccentric inner radius, m
    rolling.Rs = 10e-3                              #Shaft radius, m
    rolling.rv = 2.6e-3                             #Radius of vane tip, m
    rolling.Hc =  44.0e-3                           #Compressor cylinder height, m
    rolling.b =  4.7e-3                             #Vane thickness, m
    rolling.hv = 27.8e-3                            #Vane height, m
    rolling.L_slot = 15.2e-3                        #Length of vane slot, m
    rolling.Hshell = 100e-3                         #Height of compressor shell, m
    rolling.Rshell = 40e-3                          #Compressor shell outer diameter, m
    rolling.Nmot = 2875                             #Electric motor rotational speed, rpm
    rolling.omega = 2*pi*rolling.Nmot/60            #Frequency, rad/sec (60Hz)
    rolling.d_discharge= 8e-3                       #Discharge port diameter, m
    rolling.d_suction= 14e-3                        #Suction diameter, m
    rolling.A_discharge=pi*rolling.d_discharge**2/4
    rolling.A_suction=pi*rolling.d_suction**2/4    
    rolling.D_sp = 12e-3                   #Suction pipe diameter, m
    rolling.D_dpc = 7.5e-3                 #Diameter of processing tool for discharge cut, m
    
    # Half angle of suction port, rad
    rolling.delta_theta_s = rolling.D_sp/(2*rolling.Rc)  
    # Center angle of suction port in cylinder wrt vane, rad
    rolling.theta_sm = 22*pi/180     
    # Starting angle of suction port in cylinder wrt vane, rad
    rolling.theta_s1 = rolling.theta_sm - rolling.delta_theta_s    
    # Ending angle of suction port in cylinder wrt vane, rad
    rolling.theta_s2 = rolling.theta_sm + rolling.delta_theta_s    
    
    # Half angle of discharge port, rad
    rolling.delta_theta_d = rolling.D_dpc/(2*rolling.Rc) 
    # Center angle of discharge port in cylinder wrt vane, rad
    rolling.theta_dm_ang = 20*pi/180      
    # Center angle of discharge port in cylinder wrt vane, rad
    rolling.theta_dm = 2*pi - rolling.theta_dm_ang          
    # Starting angle of discharge port in cylinder wrt vane, rad
    rolling.theta_d1 = rolling.theta_dm - rolling.delta_theta_d  
    # Ending angle of discharge port in  cylinder wrt vane, rad 
    rolling.theta_d2 = rolling.theta_dm + rolling.delta_theta_d
    
    #Vertical height from base of tool angle to cylinder end, m
    rolling.Rp = 29e-3                     #Distance between cylinder center and end of tool angle, m
    rolling.theta_p = 59*pi/180            #Angle of processing tool discharge cut, radians
    rolling.hp = (rolling.Rp - rolling.Rc)*tan(rolling.theta_p) 
    
    # Leakage gaps
    rolling.delta_gap = 50e-6 #[m]
    rolling.delta_min = 55e-6           #Counter set clearance (at 270 degrees), m
    rolling.delta_max = 25e-6           #Set clearance (at 90 degrees), m
    rolling.delta_side = 20e-6          # Vertical clearance between roller and cylinder (one side), m
    rolling.delta_vb = 20e-6            # Vertical clearance between vane and cylinder (one side) [m]
    rolling.delta_vt = 3e-6             #Radial clearance between vane tip and roller, m
    rolling.A_ls = 692.78e-6            #Area from lower shell to upper shell, m^2
    rolling.xs = 0                      #[-]
    rolling.xc = 0                      #[-]
    rolling.x_slot = 0.5                #[-]
    rolling.x_roller = 0.8              #High concentration of oil
    
    # Flow coeff
    rolling.Cflow_roller = 0.5 #[-]
    rolling.Cflow_vb = 0.5 #[-]
    rolling.Cflow_vt = 0.5 #[-]
    rolling.Cflow_32 = 0.43 #[-]
    
    #Selection of HT model
    rolling.HT_on = True
    #These are parameters needed for the ambient heat transfer model
    rolling.A_shell = pi*rolling.Hshell*2*rolling.Rshell #[m2]
    rolling.Tamb = 25 + 273.15                    #[K] 
    
    #Parameters for the mechanical losses model (simplified)
    rolling.mech = struct()
    rolling.mech.detailed_analysis = True
    rolling.mech.Wdot_parasitic = 0.01 #Parasitic losses [kW]
    
    #Motor
    rolling.motor = Motor()
    rolling.motor.set_eta(0.8)  # Ooi, 2005
    rolling.motor.suction_fraction = 1.0  #[-]
    
    #Boundary condition
    rolling.Ref=Ref
    rolling.rho_oil = 820 #[kg/m3]
    Te_K = Te + 273.15 #[K]
    Tc_K = Tc + 273.15 
    DT_sh = 11.1 #[K]         
    Tin = Te_K + DT_sh #[K]
    DT_sc = 7 #[K]
    temp = State.State(Ref,{'T':Te_K,'Q':1}) #[K]
    pe = temp.p #[kPa]
    temp.update(dict(T=Tc_K, Q=1)) 
    pc = temp.p #[kPa]
    inletState = State.State(Ref,{'T':Tin,'P':pe})

    T2s = rolling.guess_outlet_temp(inletState,pc) #[K]
    outletState = State.State(Ref,{'T':T2s,'P':pc})
    
    #Guess mass flow rate
    rolling.Vdisp = rolling.V_disp(rolling.theta_sm,rolling.theta_dm_ang) #Displacement volume
    print('Compressor displacement (cm3/rev):',rolling.Vdisp*1e6)
    
    mdot_guess = inletState.rho*rolling.Vdisp*rolling.Nmot/60 #[kg/s]

    #Add the inlet tube
    rolling.add_tube(   Tube(key1='inlet.1',
                        key2='inlet.2',
                        L=0.03,ID=14e-3,
                        mdot=mdot_guess, 
                        State1=inletState.copy(),
                        fixed=1,
                        TubeFcn=rolling.TubeCode) 
                        )
    
    #Add the outlet tube
    rolling.add_tube( Tube(key1='outlet.1',
                            key2='outlet.2',
                            L=0.03,ID=8e-3,
                            mdot=mdot_guess, 
                            State2=outletState.copy(),
                            fixed=2,
                            TubeFcn=rolling.TubeCode) 
                            )
    
    #Add the control volumes.
    """
    'A' = Suction CV
    'B' = Compression/discharge CV
    'shell' = compressor shell
    """
    rolling.add_CV( ControlVolume(key='A',
                               initialState=inletState.copy(),
                               VdVFcn=rolling.Vs_dVs_ParkYC,
                               becomes='B') )
    rolling.add_CV( ControlVolume(key='B',
                               initialState=inletState.copy(),
                               VdVFcn=rolling.Vc_dVc_ParkYC,
                               becomes='A') )

    rolling.add_CV( ControlVolume(key='shell',
                               initialState=outletState.copy(),
                               VdVFcn=rolling.V_shell,
                               becomes='shell') )

    #Add the flow paths that link flow nodes together
    rolling.add_flow(FlowPath(key1='inlet.2',key2='A',MdotFcn=rolling.Suction))
    rolling.add_flow(FlowPath(key1='B',key2='shell',MdotFcn=rolling.Discharge))
    rolling.add_flow(FlowPath(key1='shell',key2='outlet.1',MdotFcn=rolling.Outlet))
    rolling.add_flow(FlowPath(key1='A',key2='B',MdotFcn=rolling.CompressibleGasLeakage_vb))
    rolling.add_flow(FlowPath(key1='A',key2='B',MdotFcn=rolling.CompressibleGasLeakage_vt))
    rolling.add_flow(FlowPath(key1='A',key2='B',MdotFcn=rolling.CompressibleGasLeakage_32))
    
    #Add the discharge valve parameters
    E = 1.93e11                    #Youngs Modulus, [Pa]
    h_valve = 0.2e-3               #Valve thickness, [m]
    rho_valve = 7200               #Density of spring steel, [kg/m^3] 
    C_D = 1.17                     #Drag coefficient, [-]
    d_valve_discharge = 10e-3      #Discharge valve diameter, [m]
    l_valve_discharge = 21e-3      #Total valve length, [m]
    a_valve_discharge = 16.5e-3    #Distance from anchor point to center of discharge port, [m]
    x_stopper_discharge = 3.0e-3   #Stopper height at center of discharge port, [m]

    I=(d_valve_discharge*h_valve**3)/12  #Moment of Intertia for discharge valve,[m^4]
    k_valve= (6*E*I)/(a_valve_discharge**2*(3*l_valve_discharge-a_valve_discharge))    #Valve stiffness
    m_eff=(1/3)*rho_valve*l_valve_discharge*d_valve_discharge*h_valve      #Effective mass of valve reeds, [kg]
    x_tr_discharge = 0.25*(rolling.d_discharge**2/d_valve_discharge) #Critical lift position, [m]
    
    #Define discharge valve
    rolling.discharge_valve=ValveModel(
                                        d_valve=d_valve_discharge,
                                        d_port=rolling.d_discharge,
                                        C_D= 1.14,
                                        rho_valve=rho_valve,
                                        x_stopper=x_stopper_discharge,
                                        m_eff = m_eff,
                                        k_valve = k_valve,
                                        x_tr = x_tr_discharge,
                                        key_up='B',
                                        key_down='shell'
                                        )
    #Add the discharge valve
    rolling.add_valve(rolling.discharge_valve)
    
    #Connect the callbacks for the step, endcycle, heat transfer and lump energy balance
    rolling.connect_callbacks(step_callback = rolling.step_callback,
                            endcycle_callback=rolling.endcycle_callback, 
                            heat_transfer_callback=rolling.heat_transfer_callback,
                            lumps_energy_balance_callback = rolling.lump_energy_balance_callback
                            )
    #Set debug verbosity level
    rolling.verbosity = 1
    
    t1 = timeit.default_timer()
    rolling.solver = 'RK45'
    rolling.solve(key_inlet='inlet.1',
                key_outlet='outlet.2',
                solver_method = rolling.solver,
                eps_cycle = 0.01,#0.001,
                eps_energy_balance = 0.01,
                hmin=1e-5, 
                max_number_of_steps = 100000, 
                plot_every_cycle = False
                )
    
    print('time taken', timeit.default_timer()-t1)
    
    debug_plots(rolling)
    
    del rolling.FlowStorage
    from PDSim.misc.hdf5 import HDF5Writer
    h5 = HDF5Writer()
    h5.write_to_file(rolling, HDF5file)
    
    return rolling
    
    
if __name__=='__main__':
    
    RollingPistonCompressor()
    