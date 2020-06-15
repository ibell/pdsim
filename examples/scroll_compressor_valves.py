## This file is meant for experimentation with PDSim features and may not
## work for you.  Caveat emptor!!
#
from __future__ import division, print_function
from math import pi
import timeit

from PDSim.flow.flow_models import IsentropicNozzleWrapper
from PDSim.flow.flow import FlowPath
from PDSim.core.core import struct
from PDSim.scroll.core import Scroll
from PDSim.core.containers import Tube
from PDSim.core.motor import Motor
from PDSim.flow.flow_models import ValveModel

try:
    # from PDSim.scroll.plots import plotScrollSet
    from PDSim.plot.plots import debug_plots # (Uncomment if you want to do the debug_plots)
    plotting = True
except ImportError as IE:
    print(IE)
    print('Plotting is disabled')
    plotting = False
    
from CoolProp import State

def Compressor(ScrollClass, Te = 253, Tc = 310, f = None, OneCycle = False, Ref = 'R410A', HDF5file='scroll_compressor.h5', discharge_valve=True):

    ScrollComp = ScrollClass()
    
    ScrollComp.set_scroll_geo(83e-6, 2.2, 0.005, 0.006) #Set the scroll wrap geometry
    ScrollComp.set_disc_geo('2Arc',r2 = 0)
    ScrollComp.geo.delta_flank = 10e-6
    ScrollComp.geo.delta_radial = 10e-6
    
    ScrollComp.geo.delta_suction_offset = 0.0e-3
    ScrollComp.geo.phi_ie_offset = 0.0

    # print(ScrollComp.geo); quit()
    
    ScrollComp.omega = 3000/60*2*pi
    ScrollComp.Tamb = 298.0
    
    #Temporarily set the bearing dimensions
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
    ScrollComp.mech.thrust_friction_coefficient = 0.028 #From Chen thesis
    ScrollComp.mech.orbiting_scroll_mass = 2.5
    ScrollComp.mech.L_ratio_bearings = 3
    ScrollComp.mech.mu_oil = 0.008
    
    ScrollComp.h_shell = 0.02
    ScrollComp.A_shell = 0.05
    ScrollComp.HTC = 0.01
    
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
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1 = 'sa',
                                 key2 = 's2',
                                 MdotFcn = ScrollComp.SA_S2,
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                                 )
                        )

    if not discharge_valve:
    
        ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
                                     key2 = 'dd',
                                     MdotFcn = ScrollComp.DISC_DD,
                                     MdotFcn_kwargs = dict(X_d = 0.7)
                                     )
                            )
           
        ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
                                     key2 = 'ddd',
                                     MdotFcn = ScrollComp.DISC_DD,
                                     MdotFcn_kwargs = dict(X_d = 0.7)
                                    )
                            )
    else:

        E = 1.93e11             # Youngs Modulus, [Pa]
        h_valve = 0.0006        # Valve thickness, [m]
        d_discharge = ScrollComp.geo.ra_arc1*1.9 # Port diameter [m]
        l_valve = 5*d_discharge # Total length of valve, [m]
        a_valve = l_valve/1.5   # Distance from anchor to force, [m]
        assert(a_valve < l_valve)
        rho_valve = 8000        # Density of spring steel, [kg/m^3]
        C_D = 1.17              # Drag coefficient [-]
        d_valve = d_discharge*1.5  # Valve Diameter [m]
        x_stopper = 0.006       # Stopper location [m]

        I=(d_valve*h_valve**3)/12  # Moment of Inertia for valve,[m^4]
        k_valve=(6*E*I)/(a_valve**2*(3*l_valve-a_valve))  # Valve stiffness
        m_eff=(1/3)*rho_valve*l_valve*d_valve*h_valve     # Effective mass of valve reeds
        x_tr_discharge = 0.25*(d_discharge**2/d_valve)

        # Construct the valve
        ScrollComp.discharge_valve = ValveModel(
              d_valve=d_valve,
              d_port=d_discharge,
              C_D=C_D,
              m_eff=m_eff,
              k_valve=k_valve,
              rho_valve=rho_valve, # Not used directly
              x_stopper=x_stopper,
              x_tr=x_tr_discharge,
              key_up=['ddd','dd'],
              key_down='outlet.1'
              )
        # Inform the model about the valve
        ScrollComp.add_valve(ScrollComp.discharge_valve)

        ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
                                     key2 = 'dd',
                                     MdotFcn = ScrollComp.DischargeValve,
                                     MdotFcn_kwargs = dict(X_d = 0.7)
                                     )
                            )
           
        ScrollComp.add_flow(FlowPath(key1 = 'outlet.1',
                                     key2 = 'ddd',
                                     MdotFcn = ScrollComp.DischargeValve,
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

    debug_plots(ScrollComp, family='Scroll Compressor')
    
    return ScrollComp
    
if __name__=='__main__':
    Compressor(Scroll)
