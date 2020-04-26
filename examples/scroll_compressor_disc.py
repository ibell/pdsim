"""
This file demonstrates how you can add a number of ports into the 
fixed scroll, and in this case, the ports are connecting with a discharge
plenum. 

In practice, it is probably not feasible to have so many ports, but
this example serves more as a proof-of-principle that continuous discharging
can be used to eliminate almost entirely discharge overpressure losses.

This script was based upon the computational routines developed for 
vapor injection, and then flipping the logic on its head to enable 
the opposite behavior.
"""
from __future__ import division, print_function

import sys
from math import pi
import glob

import h5py
import numpy as np
import scipy.interpolate as interp

from PDSim.flow.flow_models import IsentropicNozzleWrapper
from PDSim.flow.flow import FlowPath
from PDSim.core.core import struct
from PDSim.scroll.core import Scroll
from PDSim.core.containers import ControlVolume, Tube
from PDSim.core.motor import Motor
from PDSim.scroll.core import Port

try:
    # from PDSim.scroll.plots import plotScrollSet
    from PDSim.plot.plots import debug_plots # (Uncomment if you want to do the debug_plots)
    plotting = True
except ImportError as IE:
    print(IE)
    print('Plotting is disabled')
    plotting = False
    
from CoolProp import State
import timeit

def Compressor(ScrollClass, *, Te=273, Tc=300, OneCycle=False, Ref='R410A', HDF5file='scroll_compressor_disc.h5', Nports=4):

    ScrollComp = ScrollClass()
    #This runs if the module code is run directly 
    
    ScrollComp.set_scroll_geo(83e-6, 3.3, 0.005, 0.006) #Set the scroll wrap geometry
    ScrollComp.set_disc_geo('2Arc',r2 = 0)
    ScrollComp.geo.delta_flank = 10e-6
    ScrollComp.geo.delta_radial = 10e-6
    
    ScrollComp.geo.delta_suction_offset = 0.0e-3
    ScrollComp.geo.phi_ie_offset = 0.0
    
    ScrollComp.omega = 3000/60*2*pi
    ScrollComp.Tamb = 298.0
    
    # Temporarily set the bearing dimensions
    ScrollComp.mech = struct()
    ScrollComp.mech.D_upper_bearing = 0.04
    ScrollComp.mech.L_upper_bearing = 0.04
    ScrollComp.mech.c_upper_bearing = 20e-6
    ScrollComp.mech.D_crank_bearing = 0.04
    ScrollComp.mech.L_crank_bearing = 0.04
    ScrollComp.mech.c_crank_bearing = 20e-6
    ScrollComp.mech.D_lower_bearing = 0.025
    ScrollComp.mech.L_lower_bearing = 0.025
    ScrollComp.mech.c_lower_bearing = 20e-6
    ScrollComp.mech.thrust_ID = 0.05
    ScrollComp.mech.thrust_friction_coefficient = 0.028 # From Chen thesis
    ScrollComp.mech.orbiting_scroll_mass = 2.5
    ScrollComp.mech.L_ratio_bearings = 3
    ScrollComp.mech.mu_oil = 0.008
    
    ScrollComp.h_shell = 0.02 # kW/m^2/K
    ScrollComp.A_shell = 0.05 # m^2
    ScrollComp.HTC = 0.050 # Heat transfer coefficient between scroll and gas in chambers, in kW/m^2/K
    
    ScrollComp.motor = Motor()
    ScrollComp.motor.set_eta(0.9)
    ScrollComp.motor.suction_fraction = 1.0
    
    Tin = Te + 11.1
    temp = State.State(Ref,{'T':Te,'Q':1})
    pe = temp.p
    temp.update(dict(T=Tc, Q=1))
    pc = temp.p
    inletState = State.State(Ref,{'T':Tin,'P':pe})

    T2s = ScrollComp.guess_outlet_temp(inletState,pc)
    outletState = State.State(Ref,{'T':T2s,'P':pc})
    
    mdot_guess = inletState.rho*ScrollComp.Vdisp*ScrollComp.omega/(2*pi)
    
    ScrollComp.add_tube(Tube(key1='inlet.1',
                             key2='inlet.2',
                             L=0.3,
                             ID=0.02,
                             mdot=mdot_guess, 
                             State1=inletState.copy(),
                             fixed=1,
                             TubeFcn=ScrollComp.TubeCode))
    ScrollComp.add_tube(Tube(key1='outlet.1',
                             key2='outlet.2',
                             L=0.3,
                             ID=0.02,
                             mdot=mdot_guess, 
                             State2=outletState.copy(),
                             fixed=2,
                             TubeFcn=ScrollComp.TubeCode))
    
    ScrollComp.auto_add_CVs(inletState, outletState)
    ScrollComp.add_CV(ControlVolume(
        key='discplenum',
        initialState=outletState.copy(),
        VdVFcn=lambda theta: (1e-4, 0) # Constant volume, outputs of function are volume, dV/dtheta
    ))
    
    ScrollComp.auto_add_leakage(flankFunc = ScrollComp.FlankLeakage, 
                                radialFunc = ScrollComp.RadialLeakage)
    
    FP = FlowPath(key1='inlet.2', 
                  key2='sa', 
                  MdotFcn=IsentropicNozzleWrapper(),
                  )
    FP.A = pi*0.01**2/4
    ScrollComp.add_flow(FP)
    FP = FlowPath(key1='outlet.1', 
                  key2='discplenum', 
                  MdotFcn=IsentropicNozzleWrapper(),
                  )
    FP.A = pi*0.01**2/4
    ScrollComp.add_flow(FP)
    
    ScrollComp.add_flow(FlowPath(key1='sa', 
                                 key2='s1',
                                 MdotFcn=ScrollComp.SA_S1,
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1 = 'sa',
                                 key2 = 's2',
                                 MdotFcn = ScrollComp.SA_S2,
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                                 )
                        )
    
    ScrollComp.add_flow(FlowPath(key1 = 'discplenum',
                                 key2 = 'dd',
                                 MdotFcn = ScrollComp.DISC_DD,
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                                 )
                        )
       
    ScrollComp.add_flow(FlowPath(key1 = 'discplenum',
                                 key2 = 'ddd',
                                 MdotFcn = ScrollComp.DISC_DD,
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                                 )
                        )
#     ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
#                                  key2 = 'd1',
#                                  MdotFcn = ScrollComp.DISC_D1,
#                                  MdotFcn_kwargs = dict(X_d = 0.7)
#                                  )
#                         )
#     
#     FP = FlowPath(key1='outlet.1', 
#                   key2='dd', 
#                   MdotFcn=IsentropicNozzleWrapper(),
#                   )
#     FP.A = pi*0.006**2/4
#     ScrollComp.add_flow(FP)
#       
#     FP = FlowPath(key1='outlet.1', 
#                   key2='ddd', 
#                   MdotFcn=IsentropicNozzleWrapper(),
#                   )
#     FP.A = pi*0.006**2/4
#     ScrollComp.add_flow(FP)

    # ***************
    # Discharge Ports
    # ***************

    sim = ScrollComp
    for involute in ['i', 'o']:

        # Define the parameters for the ports
        if involute == 'i':
            phi_list = np.linspace(sim.geo.phi_fie-2*np.pi-0.1, sim.geo.phi_fis+np.pi/2, Nports)
        elif involute == 'o':
            phi_list = np.linspace(sim.geo.phi_foe-2*np.pi-0.1, sim.geo.phi_fos+np.pi/2, Nports)
        else:
            raise ValueError(involute)
        involute_list = [involute]*len(phi_list)
        D_list = [sim.geo.t]*len(phi_list) # All ports are thickness of scroll in diameter
        offset_list = [d*0.5 for d in D_list] # Distance of translation of center of port normal to wall
        X_d_list = [1.0]*len(phi_list) # Flow coefficient for forward flow
        X_d_backflow_list = [0.0]*len(phi_list) # Backflow coefficient

        # Add the ports to the model
        sim.fixed_scroll_ports = []
        for i in range(Nports):
            p = Port()
            p.phi = phi_list[i]
            p.involute = involute_list[i]
            p.offset = offset_list[i]
            p.D = D_list[i]
            p.parent = 'discplenum'
            p.X_d = X_d_list[i]
            p.X_d_backflow = X_d_backflow_list[i]
            sim.fixed_scroll_ports.append(p)

        # Calculate the areas between each port and every control volume
        # 
        # Can turn on plotting (plot=True) to see the flow areas.
        sim.calculate_port_areas(plot=False)
            
        # Iterate over the ports, connecting them to discharge plenum
        for port in sim.fixed_scroll_ports:

            # Iterate over the chambers each port is connected to at least
            # at some point
            for partner in port.area_dict:
            
                # Create a spline interpolator object for the area between port and the partner chamber
                A_interpolator = interp.splrep(port.theta, port.area_dict[partner], k=2, s=0)
                
                # Add the flow path connecting the working chamber and discharge plenum
                #
                # It is a simplified nozzle model, allowing for flow out the nozzle, 
                # according to the flow area at the given crank angle, but no backflow.
                # This represents an idealized case.
                sim.add_flow(FlowPath(key1='discplenum',
                                      key2=partner,
                                      MdotFcn=sim.INTERPOLATING_NOZZLE_FLOW,
                                      MdotFcn_kwargs=dict(X_d=port.X_d,
                                                          X_d_backflow=port.X_d_backflow,
                                                          upstream_key=partner,
                                                          A_interpolator=A_interpolator
                                                          )
                                     )
                             )
    
    ScrollComp.add_flow(FlowPath(key1='d1',
                                 key2='dd',
                                 MdotFcn=ScrollComp.D_to_DD))
    ScrollComp.add_flow(FlowPath(key1='d2',
                                 key2='dd',
                                 MdotFcn=ScrollComp.D_to_DD))
    
    #Connect the callbacks for the step, endcycle, heat transfer and lump energy balance
    ScrollComp.connect_callbacks(step_callback = ScrollComp.step_callback,
                                 endcycle_callback = ScrollComp.endcycle_callback,
                                 heat_transfer_callback = ScrollComp.heat_transfer_callback,
                                 lumps_energy_balance_callback = ScrollComp.lump_energy_balance_callback
                                 )
    
    
    t1 = timeit.default_timer()
    ScrollComp.RK45_eps = 1e-6
    ScrollComp.eps_cycle = 3e-3
    try:
        ScrollComp.solve(key_inlet='inlet.1',
                         key_outlet='outlet.2',
                         solver_method='RK45',
                         OneCycle = OneCycle,
                         plot_every_cycle= False,
                         #hmin = 1e-3
                         eps_cycle = 3e-3
                         )
    except BaseException as E:
        print(E)
        raise

    print('time taken', timeit.default_timer()-t1)
    
    del ScrollComp.FlowStorage
    from PDSim.misc.hdf5 import HDF5Writer
    h5 = HDF5Writer()
    h5.write_to_file(ScrollComp, HDF5file)

    if '--plot' in sys.argv:
        debug_plots(ScrollComp, family='Scroll Compressor')
    
    return ScrollComp
    
if __name__=='__main__':

    for Nports in [0, 2, 3, 4, 5, 6]:
        Compressor(Scroll, HDF5file='scroll_compressor_'+str(Nports)+'disc.h5', Nports=Nports)
    
    for fname in glob.glob('scroll_compressor_*disc.h5'):
        with h5py.File(fname) as f:
            print(fname, f['eta_oi'][()], f['eta_a'][()])