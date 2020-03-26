## This file is meant for experimentation with PDSim features and may not
## work for you.  Caveat emptor!!
#
from __future__ import division, print_function

import sys, os
from math import pi

from PDSim.flow.flow_models import IsentropicNozzleWrapper
from PDSim.flow.flow import FlowPath
from PDSim.scroll import scroll_geo
from PDSim.core.core import struct
from PDSim.scroll.core import Scroll
from PDSim.core.containers import ControlVolume, Tube
from PDSim.core.motor import Motor
try:
    from PDSim.scroll.plots import plotScrollSet
    from PDSim.plot.plots import debug_plots # (Uncomment if you want to do the debug_plots)
    plotting = True
except ImportError as IE:
    print(IE)
    print('Plotting is disabled')
    plotting = False

from CoolProp import State
from CoolProp import CoolProp as CP

import numpy as np
import time

def Compressor(Te = 0,DTsh = 11.1,Tc = 20, Tamb = 25, Nmot = 3600, f = None, OneCycle = False, Ref = 'R134a', HDF5file='scroll_compressor.h5'):

    ScrollComp=Scroll()

    #Set the scroll wrap geometry
    ScrollComp.set_scroll_geo(82.7e-6,2.5, 0.004,0.005,phi_i0=0.0,phi_os=0.3, phi_is = 3.142)
    ScrollComp.set_disc_geo('2Arc',r2 = 0.0010)

    ScrollComp.geo.delta_suction_offset = 0.0e-3
    ScrollComp.geo.phi_ie_offset = 0.0

    ScrollComp.omega = Nmot/60*2*pi
    ScrollComp.slip_ratio = 0.01
    ScrollComp.Tamb = Tamb + 273.15 #[K]

    #Temporarily set the bearing dimensions
    ScrollComp.mech = struct()
    ScrollComp.mech.detailed_analysis = True
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
    ScrollComp.mech.thrust_friction_coefficient = 0.013#0.028 #From Chen thesis
    ScrollComp.mech.orbiting_scroll_mass = 2.5
    ScrollComp.mech.L_ratio_bearings = 3
    ScrollComp.mech.mu_oil = 0.005
    ScrollComp.mech.rho_oil = 930

    ScrollComp.mech.oldham_ring_radius = 0.06 #m
    ScrollComp.mech.oldham_mass = 0.1 #kg
    ScrollComp.mech.oldham_thickness = 0.008 #m
    ScrollComp.mech.oldham_key_height = 0.006 #m
    ScrollComp.mech.oldham_key_width = 0.006 #m
    ScrollComp.mech.oldham_key_friction_coefficient = 0.01 #-
    ScrollComp.mech.oldham_rotation_beta = 0 #rad

    ScrollComp.h_shell = 0.02
    ScrollComp.A_shell = 450.75e-3*((246.126e-3)**2*pi/4)
    ScrollComp.HTC = 0.0
    ScrollComp.HT_corr = 'Pereira-Deschamps' #'Jang-Jeong'

    # Temperature Lumps
    ScrollComp.OEB_type = 'multi-lump' #'single-lump'
    ScrollComp.OEB_solver = 'MDNR'
    ScrollComp.Rshell_oil = 190 #K/kW  from Chen (2000) - PhD thesis

    # Define motor efficiency
    ScrollComp.motor.set_eta(0.95)
    ScrollComp.motor.suction_fraction = 0.5

    # Operating conditions
    Tevap = Te + 273.15 #K
    Tcond = Tc + 273.15 #K

    Tin = Tevap + DTsh
    DT_sc = 7 #K
    temp = State.State(Ref,{'T':Tevap,'Q':1})
    pe = temp.p
    temp.update(dict(T=Tcond, Q=1))
    pc = temp.p
    inletState = State.State(Ref,{'T':Tin,'P':pe})
    ScrollComp.pressure_ratio = pc/pe
    print('Pressure Ratio:',pc/pe, pc,pe )
    T2s = ScrollComp.guess_outlet_temp(inletState,pc)
    outletState = State.State(Ref,{'T':T2s,'P':pc})
    print('Vdisp:', ScrollComp.Vdisp*1e6,'cm3/rev')
    mdot_guess = inletState.rho*ScrollComp.Vdisp*ScrollComp.omega/(2*pi)

    #Leakage flow
    ScrollComp.geo.delta_flank = -9.615e-7*(ScrollComp.pressure_ratio-1.67) + 15e-6  #10e-6
    ScrollComp.geo.delta_radial = 1.1e-6*(ScrollComp.pressure_ratio-1.67) + 1e-6  #10e-6


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

    ScrollComp.auto_add_leakage(flankFunc = ScrollComp.FlankLeakage,
                                radialFunc = ScrollComp.RadialLeakage)

    FP = FlowPath(key1='inlet.2',
                  key2='sa',
                  MdotFcn=IsentropicNozzleWrapper(),
                  )
    FP.A = pi*0.01**2/4
    ScrollComp.add_flow(FP)

    ScrollComp.add_flow(FlowPath(key1='sa',
                                 key2='s1',
                                 MdotFcn=ScrollComp.SA_S1,
                                 MdotFcn_kwargs = dict(X_d = 0.8)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1 = 'sa',
                                 key2 = 's2',
                                 MdotFcn = ScrollComp.SA_S2,
                                 MdotFcn_kwargs = dict(X_d = 0.8)
                                 )
                        )

    ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
                                 key2 = 'dd',
                                 MdotFcn = ScrollComp.DISC_DD,
                                 MdotFcn_kwargs = dict(X_d = 0.8)
                                 )
                        )

    ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
                                 key2 = 'ddd',
                                 MdotFcn = ScrollComp.DISC_DD,
                                 MdotFcn_kwargs = dict(X_d = 0.8)
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

    from time import clock
    t1=clock()
    ScrollComp.RK45_eps = 1e-6
    ScrollComp.eps_cycle = 3e-3
    ScrollComp.verbosity = 10

    ScrollComp.precond_solve(key_inlet='inlet.1',
                                key_outlet='outlet.2',
                                solver_method='RK45',
                                OneCycle = OneCycle,
                                plot_every_cycle= False,
                                x0 = [330,330,350], #Guesses [Td,Tlump[0],Tlump[1]]
                                #hmin = 1e-3
                                eps_cycle = 3e-3,
                                eps_energy_balance = 0.1 #relaxed multi-lump convergence
                                )

    print('time taken', clock()-t1)

    #debug_plots(ScrollComp)

    del ScrollComp.FlowStorage
    from PDSim.misc.hdf5 import HDF5Writer
    h5 = HDF5Writer()
    h5.write_to_file(ScrollComp, HDF5file)

    return ScrollComp


if __name__=='__main__':

    # Example
    scroll = Compressor(Te = 0,
                        DTsh = 11.1,
                        Tc = 20,
                        Tamb = 25,
                        Nmot = 3600,
                        f = None,
                        OneCycle = False,
                        Ref = 'R134a',
                        HDF5file='ScrollCompressor_MultiLump.h5')


