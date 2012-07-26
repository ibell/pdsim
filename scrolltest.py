
# If being run from the folder that contains the PDSim source tree, 
# remove the current location from the python path and use the 
# site-packages version of PDSim
import sys, os

try:
    # if you can import listmath
    import PDSim.misc._listmath
except ImportError:
    #If you can't import listmath, remove the PDSim folder
    current_path = os.path.abspath(os.curdir)
    if current_path in sys.path:
        i = sys.path.index(current_path)
        sys.path.pop(i)
    else:
        raise ImportError('_listmath was not found and PDSim folder not found on the path') 

from math import pi

from PDSim.flow.flow import FlowPath
from PDSim.scroll import scroll_geo
from PDSim.scroll.core import Scroll
from PDSim.core.containers import ControlVolume
from PDSim.core.core import Tube
from PDSim.plot.plots import debug_plots
from PDSim.scroll.plots import plotScrollSet

from CoolProp import State
from CoolProp import CoolProp as CP

import numpy as np

from matplotlib import pyplot as plt
import time

def Compressor():
    
    ScrollComp=Scroll()
    #This runs if the module code is run directly
    ScrollComp.set_scroll_geo(104.8e-6,2.2,0.004,0.005) #Set the scroll wrap geometry
    ScrollComp.set_disc_geo('2Arc',r2='PMP')
    ScrollComp.geo.delta_flank=15e-6
    ScrollComp.geo.delta_radial=15e-6
    ScrollComp.omega=3500/60*2*pi
    ScrollComp.Tamb = 298.0
    
    scroll_geo.plot_HT_angles(pi, ScrollComp.geo, ['s1','s2','c1.1','c1.2'], 'i')
    
    Injection = False
        
    Ref='R410A'
    #State.debug(10)
    State.set_1phase_LUT_params(Ref,10,10,250,500,200,3000)
    State.LUT(True)
    
    inletState = State.State(Ref,{'T':300,'P':310})
    outletState = State.State(Ref,{'T':400,'P':1200})
    
    mdot_guess = inletState.rho*ScrollComp.Vdisp*ScrollComp.omega/(2*pi)
    
    ScrollComp.add_tube(Tube(key1='inlet.1',key2='inlet.2',L=0.3,ID=0.02,
                             mdot=mdot_guess, State1=inletState.copy(),
                             fixed=1,TubeFcn=ScrollComp.TubeCode))
    ScrollComp.add_tube(Tube(key1='outlet.1',key2='outlet.2',L=0.3,ID=0.02,
                             mdot=mdot_guess, State2=outletState.copy(),
                             fixed=2,TubeFcn=ScrollComp.TubeCode))
    
    if Injection:
        theta = 0
#        #Plot the injection ports (symmetric)
#        plotScrollSet(theta, ScrollComp.geo)
#        ax = plt.gca()
        phi = ScrollComp.geo.phi_oe-pi-2*pi+0.01
#        #Involute angle along the outer involute of the scroll wrap
#        x,y = scroll_geo.coords_inv(phi,ScrollComp.geo,theta,'fo')
#        nx,ny = scroll_geo.coords_norm(phi,ScrollComp.geo,theta,'fo')
#        rport = 0.002
#        xc,yc = x-nx*rport,y-ny*rport
#        ax.plot(xc, yc, '.')
#        t = np.linspace(0,2*pi,100)
#        ax.plot(xc + rport*np.cos(t),yc+rport*np.sin(t),'k')
#        plt.show()
#        
#        for theta in np.linspace(0,2*pi,15):
#            print theta, ScrollComp._get_injection_CVkey(phi + pi, theta, inner_outer='i')
#        for theta in np.linspace(0,2*pi,15):
#            print theta, ScrollComp._get_injection_CVkey(phi, theta, inner_outer='o')
        
    if Injection:
        #Tube is a meter long with ID of 0.01 m
        f = 0.6
        p = f*outletState.p+(1-f)*inletState.p
        Tsat = CP.Props('T', 'P', p, 'Q', 1.0, Ref)
        T = Tsat + 3.0
        rho = CP.Props('D', 'T', T, 'P',p, Ref)
        injState1 = State.State(Ref, dict(T=T, D=rho))
        V_tube = 1.0*pi*0.01**2/4.0
        ScrollComp.add_CV(ControlVolume(key ='injCV.1',
                                        VdVFcn = ScrollComp.V_injection,
                                        VdVFcn_kwargs = dict(V_tube=V_tube),
                                        initialState = injState1,
                                        becomes = 'injCV.1'
                                        )
                          )
    
    ScrollComp.auto_add_CVs(inletState, outletState)
    ScrollComp.auto_add_leakage(flankFunc = ScrollComp.FlankLeakage, 
                                radialFunc = ScrollComp.RadialLeakage)
    
    print "still missing radial leakage terms"
    
    ScrollComp.add_flow(FlowPath(key1='inlet.2', key2='sa', MdotFcn=ScrollComp.Inlet_sa,
                                 MdotFcn_kwargs = dict(X_d = 0.4)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1='sa', key2='s1', MdotFcn=ScrollComp.SA_S,
                                 MdotFcn_kwargs = dict(X_d = 0.3)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1='sa', key2='s2', MdotFcn=ScrollComp.SA_S,
                                 MdotFcn_kwargs = dict(X_d = 0.3)
                                 )
                        )
    
    if Injection:
        #Injection flow paths
        ScrollComp.add_flow(FlowPath(key1= 'c1.1', 
                                     key2 = 'injCV.1', 
                                     MdotFcn=ScrollComp.Injection_to_Comp,
                                     MdotFcn_kwargs = dict(phi = phi + pi,
                                                           inner_outer = 'i')
                                    )
                            )
        ScrollComp.add_flow(FlowPath(key1 = 'c2.1', 
                                     key2 = 'injCV.1', 
                                     MdotFcn=ScrollComp.Injection_to_Comp,
                                     MdotFcn_kwargs = dict(phi = phi,
                                                           inner_outer = 'o')
                                    )
                            )
        
        ScrollComp.add_tube(Tube(key1='injection.1',key2='injection.2',
                                 L=0.3,ID=0.02,
                                 mdot=mdot_guess, 
                                 State1=ScrollComp.CVs['injCV.1'].State.copy(),
                                 fixed=1,
                                 TubeFcn=ScrollComp.TubeCode
                                 )
                            )
        ScrollComp.add_flow(FlowPath(key1='injection.2',
                                     key2='injCV.1',
                                     MdotFcn=ScrollComp.InjectionTubeFM))
    
    ScrollComp.add_flow(FlowPath(key1='outlet.1',key2='ddd',MdotFcn=ScrollComp.Discharge))
    ScrollComp.add_flow(FlowPath(key1='outlet.1',key2='dd',MdotFcn=ScrollComp.Discharge))
    ScrollComp.add_flow(FlowPath(key1='d1',key2='dd',MdotFcn=ScrollComp.D_to_DD))
    ScrollComp.add_flow(FlowPath(key1='d2',key2='dd',MdotFcn=ScrollComp.D_to_DD))

    from time import clock
    t1=clock()
    ScrollComp.RK45_eps = 1e-8
    ScrollComp.EulerN = 20000
    ScrollComp.HeunN = 6000
    ScrollComp.solve(key_inlet='inlet.1',key_outlet='outlet.2',
                 step_callback=ScrollComp.step_callback, 
                 endcycle_callback=ScrollComp.endcycle_callback,
                 heat_transfer_callback=ScrollComp.heat_transfer_callback,
                 lump_energy_balance_callback=ScrollComp.lump_energy_balance_callback,
                 solver_method='RK45',
                 hmin=2*pi/(100000000),
                 UseNR = False, #Use Newton-Raphson ND solver to determine the initial state
                 OneCycle = False
                 )
    
    print 'time taken',clock()-t1
    
    if Injection:
        print ScrollComp.FlowsProcessed.mean_mdot
    debug_plots(ScrollComp)

    
if __name__=='__main__':
        
    profile=False
    if profile==True:
        import line_profiler as LP
        profiler=LP.LineProfiler(Scroll.cycle_RK45)
        profiler.run("Compressor()")
        profiler.print_stats()
    else:
        Compressor()

