## This file is meant for experimentation with PDSim features and may not
## work for you.  Caveat emptor!!
#
from __future__ import division

# If being run from the folder that contains the PDSim source tree, 
# remove the current location from the python path and use the 
# site-packages version of PDSim
import sys, os
from math import pi

# If the following line is uncommented, python will try to use a local version
# of PDSim.  This is handy for debugging purposes.  Generally you want this line 
# commented out
#sys.path.insert(0, os.path.abspath('..'))

from PDSim.flow.flow import FlowPath
from PDSim.scroll import scroll_geo
from PDSim.scroll.core import Scroll
from PDSim.core.containers import ControlVolume
from PDSim.core.core import Tube
from PDSim.plot.plots import debug_plots
from PDSim.scroll.plots import plotScrollSet
from PDSim.core.motor import Motor

from CoolProp import State
from CoolProp import CoolProp as CP

import numpy as np

from matplotlib import pyplot as plt
import time

Injection = False
check_valve = False

def Compressor(f = None):
    global Injection
    ScrollComp=Scroll()
    #This runs if the module code is run directly
    ScrollComp.set_scroll_geo(104.8e-6, 1.61, 0.004, 0.005) #Set the scroll wrap geometry
    ScrollComp.set_disc_geo('2Arc',r2='PMP')
    ScrollComp.geo.delta_flank = 1.5e-6
    ScrollComp.geo.delta_radial = 1.5e-6
    ScrollComp.omega = 3600/60*2*pi
    ScrollComp.Tamb = 298.0
    ScrollComp.eta_motor = 0.9
    
    #Temporarily set the bearing dimensions
    ScrollComp.D_upper_bearing = 0.025
    ScrollComp.L_upper_bearing = 0.025
    ScrollComp.c_upper_bearing = 20e-6
    ScrollComp.D_crank_bearing = 0.025
    ScrollComp.L_crank_bearing = 0.025
    ScrollComp.c_crank_bearing = 20e-6
    ScrollComp.D_lower_bearing = 0.025
    ScrollComp.L_lower_bearing = 0.025
    ScrollComp.c_lower_bearing = 20e-6
    ScrollComp.thrust_friction_coefficient = 0.028 #From Chen thesis
    ScrollComp.mu_oil = 0.0086
    ScrollComp.L_ratio_bearings = 3.0
    
    ScrollComp.h_shell = 10
    ScrollComp.A_shell = 0.05
    ScrollComp.HTC = 0.0
    
    ScrollComp.geo.delta_suction_offset = 0.0e-3
    ScrollComp.geo.phi_ie_offset = pi
    
    ScrollComp.motor = Motor()
    ScrollComp.motor.set_eta(0.9)
    
    ScrollComp.orbiting_scroll_mass = 3.5
    
#    print ScrollComp.V_s1(0)[0]
#    print ScrollComp.V_sa(2*pi)[0]-ScrollComp.V_sa(0)[0]
#    
#    plotScrollSet(2*pi,ScrollComp.geo,shaveOn=False,offsetScroll=ScrollComp.geo.phi_ie_offset>0,show=True)
#    theta = np.linspace(0,2*pi,101)
#    _Vs1 = []
#    _Vs2 = []
#    _Vs3 = [] 
#    for th in theta:
#        _Vs1.append(scroll_geo.S1(th,ScrollComp.geo)[0])
#    ScrollComp.geo.phi_ie_offset = pi
#    for th in theta:
#        _Vs2.append(scroll_geo.S1(th,ScrollComp.geo)[0])
#    _theta = np.linspace(0,pi,101)
##    for th in _theta:
##        _Vs3.append(scroll_geo.S1(th,ScrollComp.geo,poly=True)[2])
#
#    import pylab
#    if len(_Vs3)>0:
#        pylab.plot(theta,_Vs1,theta,_Vs2,_theta,_Vs3)
#    else:
#        pylab.plot(theta,_Vs1,theta,_Vs2)
#    pylab.show()
#    return
#        
    
#    ScrollComp.geo.phi_ie_offset = pi
#    for th in np.linspace(0,2*pi,11):
#        plotScrollSet(th,ScrollComp.geo,shaveOn=False,offsetScroll=ScrollComp.geo.phi_ie_offset>0,show=True)
#        print scroll_geo.S1(th,ScrollComp.geo,poly=True)
#    return
    
#    import pylab
#    pylab.plot(fx,fy,'-',fxp,fyp,'s',mfc='none')
#    pylab.show()
#    return
#    radial_pairs = scroll_geo.radial_leakage_pairs(ScrollComp.geo)
#    th = 0.55*2*pi
#    
#    plotScrollSet(th,ScrollComp.geo)
#    ax = plt.gca()
#    for key1, key2 in radial_pairs:
#        try:
#            phi_min, phi_max = scroll_geo.radial_leakage_angles(th, ScrollComp.geo,key1,key2)
#            phi = np.linspace(phi_min, phi_max, 50)
#            print key1, key2, phi_min, phi_max
#            
#            x,y = scroll_geo.coords_inv(phi, ScrollComp.geo, th, flag="fi")
#            ax.plot(x,y)
#            ax.plot(x[0],y[0],'o')
#            ax.plot(x[-1],y[-1],'o')
#            x,y = scroll_geo.coords_inv(phi[len(phi)//2], ScrollComp.geo, th, flag="fi")
#            ax.text(x,y,key1+'-'+key2,ha='center',va='center')
#            
#        except KeyError:
#            print 'no match for',key1,key2
#            pass
#    plt.show() 
    
    if f is None:
        Injection = False
        
    Ref='R404A'
    
    Te = -10 + 273.15
    Tc =  43 + 273.15
    Tin = Tc + 20
    DT_sc = 7
    pe = CP.Props('P','T',Te,'Q',1.0,Ref)
    pc = CP.Props('P','T',Tc,'Q',1.0,Ref)
    inletState = State.State(Ref,{'T':Tin,'P':pe})
    T2s = ScrollComp.guess_outlet_temp(inletState,pc)
    outletState = State.State(Ref,{'T':T2s,'P':pc})

#    inletState = State.State(Ref,{'T':300.0,'P':300.0})
#    p_outlet = inletState.p*4.0
#    T2s = ScrollComp.isentropic_outlet_temp(inletState,p_outlet)
#    outletState = State.State(Ref,{'T':T2s,'P':p_outlet})    
    
    mdot_guess = inletState.rho*ScrollComp.Vdisp*ScrollComp.omega/(2*pi)
    
    ScrollComp.add_tube(Tube(key1='inlet.1',key2='inlet.2',L=0.3,ID=0.02,
                             mdot=mdot_guess, State1=inletState.copy(),
                             fixed=1,TubeFcn=ScrollComp.TubeCode))
    ScrollComp.add_tube(Tube(key1='outlet.1',key2='outlet.2',L=0.3,ID=0.02,
                             mdot=mdot_guess, State2=outletState.copy(),
                             fixed=2,TubeFcn=ScrollComp.TubeCode))
             
    if Injection:
        phi = ScrollComp.geo.phi_oe-pi-2*pi+0.01
        #Tube is a meter long with ID of 0.01 m
        p = f*outletState.p+(1-f)*inletState.p
        Tsat = CP.Props('T', 'P', p, 'Q', 1.0, Ref)
        T = Tsat + 3.0
        rho = CP.Props('D', 'T', T, 'P',p, Ref)
        injState1 = State.State(Ref, dict(T=T, D=rho))
        V_tube = 1.0*pi*0.01**2/4.0
        ScrollComp.add_CV(ControlVolume(key ='injCV.1',
                                        VdVFcn = ScrollComp.V_injection,
                                        VdVFcn_kwargs = dict(V_tube = V_tube),
                                        initialState = injState1
                                        )
                          )
    
    ScrollComp.auto_add_CVs(inletState, outletState)
    ScrollComp.auto_add_leakage(flankFunc = ScrollComp.FlankLeakage, 
                                radialFunc = ScrollComp.RadialLeakage)
    
    ScrollComp.add_flow(FlowPath(key1='inlet.2', 
                                 key2='sa', 
                                 MdotFcn=ScrollComp.IsentropicNozzleFM,
                                 MdotFcn_kwargs = dict(A = pi*0.02**2/4)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1='sa', 
                                 key2='s1',
                                 MdotFcn=ScrollComp.SA_S1,
                                 MdotFcn_kwargs = dict(X_d = 0.3)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1='sa',
                                 key2='s2',
                                 MdotFcn=ScrollComp.SA_S2,
                                 MdotFcn_kwargs = dict(X_d = 0.3)
                                 )
                        )
    
    if Injection:
        #Injection flow paths
        ScrollComp.add_flow(FlowPath(key1= 'c1.1', 
                                     key2 = 'injCV.1', 
                                     MdotFcn=ScrollComp.Injection_to_Comp,
                                     MdotFcn_kwargs = dict(phi = phi + pi,
                                                           inner_outer = 'i',
                                                           check_valve = check_valve)
                                    )
                            )
        ScrollComp.add_flow(FlowPath(key1 = 'c2.1', 
                                     key2 = 'injCV.1', 
                                     MdotFcn=ScrollComp.Injection_to_Comp,
                                     MdotFcn_kwargs = dict(phi = phi,
                                                           inner_outer = 'o',
                                                           check_valve = check_valve)
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
                                     MdotFcn=ScrollComp.IsentropicNozzleFM,
                                     MdotFcn_kwargs = dict(A=pi*0.02**2/4)
                                     )
                            )
    
    ScrollComp.add_flow(FlowPath(key1='outlet.1',key2='ddd',
                                 MdotFcn = ScrollComp.IsentropicNozzleFM,
                                 MdotFcn_kwargs = dict(A=pi*0.01**2/4)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1='outlet.1',key2='dd',
                                 MdotFcn = ScrollComp.IsentropicNozzleFM,
                                 MdotFcn_kwargs = dict(A=pi*0.01**2/4)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1='d1',key2='dd',MdotFcn=ScrollComp.D_to_DD))
    ScrollComp.add_flow(FlowPath(key1='d2',key2='dd',MdotFcn=ScrollComp.D_to_DD))

    from time import clock
    t1=clock()
    ScrollComp.RK45_eps = 1e-7
    ScrollComp.EulerN = 20000
    ScrollComp.HeunN = 6000
    try:
        ScrollComp.precond_solve(key_inlet='inlet.1',key_outlet='outlet.2',
                     step_callback=ScrollComp.step_callback, 
                     endcycle_callback=ScrollComp.endcycle_callback,
                     heat_transfer_callback=ScrollComp.heat_transfer_callback,
                     lump_energy_balance_callback=ScrollComp.lump_energy_balance_callback,
                     solver_method='RK45',
                     UseNR = False, #Use Newton-Raphson ND solver to determine the initial state
                     OneCycle = True,
                     plot_every_cycle= False
                     )
    except:
        
        #debug_plots(ScrollComp)
        raise

    print 'time taken',clock()-t1
    
    if Injection:
        print ScrollComp.FlowsProcessed.mean_mdot
        ha = CP.Props('H','T',Tc-DT_sc,'P',pc,'R404A')
        hb = CP.Props('H','T',Tsat,'Q',0.0,'R404A')
        hc = CP.Props('H','T',Tsat,'Q',1.0,'R404A')
        
        ScrollComp.injection_massflow_ratio = (ha-hb)/(hc-ha)
        print 'enthalpies',ha,hb,hc,'x',ScrollComp.injection_massflow_ratio
    
    #return 
    debug_plots(ScrollComp, plot_names=['Pressure v. crank angle'])
    
    ScrollComp.calculate_force_terms(orbiting_back_pressure = pe)
    
    import pylab
    V = ScrollComp.geo.ro*ScrollComp.omega
    mu = 0.08
    print 'meanFz',ScrollComp.forces.mean_Fz
    Wdot_loss_thrust = ScrollComp.forces.mean_Fz*V*mu
    pylab.plot(ScrollComp.t,ScrollComp.forces.Fz.T,
               ScrollComp.t,ScrollComp.forces.summed_Fz,'-',
               ScrollComp.t,ScrollComp.forces.mean_Fz*(1+0*ScrollComp.t),'--') #kN
    pylab.show()
    print 'Thrust bearing loss is',Wdot_loss_thrust*1000,'W'
    
    print 'meanFr',ScrollComp.forces.mean_Fr, type(ScrollComp.forces.mean_Fr)
    pylab.plot(ScrollComp.t,ScrollComp.forces.Fr.T, #.T for transpose of the matrix
               ScrollComp.t,ScrollComp.forces.mean_Fr*(1+0*ScrollComp.t),'--') #kN
    pylab.show()
    
    from PDSim.core.bearings import journal_bearing
    JB = journal_bearing(r_b = 0.02245/2,
                    L = 0.027,
                    omega = ScrollComp.omega,
                    W = ScrollComp.forces.mean_Fr*1000,
                    design = 'friction',
                    eta_0 = 0.008
                    )
    
    print 'Crank pin journal loss is',JB['Wdot_loss'],'W'
    print 'Crank pin gap width is',JB['c']*1e6,'um'
    print 'Bearing number',JB['S']
    
    JB = journal_bearing(r_b = 0.01,
                         L = 0.02,
                         omega = ScrollComp.omega,
                         W = ScrollComp.forces.mean_Fr*1000,
                         design = 'friction',
                         eta_0 = 0.008
                         )
    print 'Upper journal loss is',JB['Wdot_loss'],'W'
    print 'Upper journal gap width is',JB['c']*1e6,'um'
        
    return ScrollComp
    
if __name__=='__main__':
    
    profile=False
    if profile==True:
        import line_profiler as LP
        profiler=LP.LineProfiler(Scroll.cycle_RK45)
        profiler.run("Compressor()")
        profiler.print_stats()
    else:
        if Injection:
            FP = open('results-checkvalve-'+str(check_valve)+'.csv','w')
            FP.write('f,mdot_injection,mdot_suction,x_sim,x_cycle\n')
            for f in np.linspace(0.2,0.8,7):
                S = Compressor(f)
                x_sim = str(S.FlowsProcessed.mean_mdot['injection.1']/S.FlowsProcessed.mean_mdot['inlet.1']) 
                x_cycle = str(S.injection_massflow_ratio)
                FP.write(str(f)+','+str(S.FlowsProcessed.mean_mdot['injection.1'])+','+str(S.FlowsProcessed.mean_mdot['inlet.1'])+','+x_sim+','+x_cycle+'\n')
            FP.close()       
        else:
            Compressor()

