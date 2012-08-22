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

from CoolProp import State
from CoolProp import CoolProp as CP

import numpy as np

from matplotlib import pyplot as plt
import time

Injection = False

def Compressor(f = None):
    global Injection
    ScrollComp=Scroll()
    #This runs if the module code is run directly
    ScrollComp.set_scroll_geo(104.8e-6, 1.61, 0.004, 0.005) #Set the scroll wrap geometry
    ScrollComp.set_disc_geo('2Arc',r2='PMP')
    ScrollComp.geo.delta_flank = 15e-6
    ScrollComp.geo.delta_radial = 15e-6
    ScrollComp.omega = 3600/60*2*pi
    ScrollComp.Tamb = 298.0
    
    M_O,cx,cy,fx,fy,fz,tt,fyp,fxp=[],[],[],[],[],[],[],[],[]
        
    for th in np.linspace(0,2*pi,200):
        tt.append(th)
        F = scroll_geo.S1_forces(th,geo=ScrollComp.geo,poly=True)
        fxp.append(F['fxp_poly'])
        fyp.append(F['fyp_poly'])
        fx.append(F['fx_p'])
        fy.append(F['fy_p'])
        fz.append(F['fz_p'])
        cx.append(F['cx'])
        cy.append(F['cy'])
        M_O.append(F['M_O'])
            
    for th in np.linspace(0,scroll_geo.theta_d(ScrollComp.geo)-0.001,100):
        tt.append(th+2*pi)
        F = scroll_geo.C1_forces(th,alpha=1,geo=ScrollComp.geo,poly=True)
        fxp.append(F['fxp_poly'])
        fyp.append(F['fyp_poly'])
        fx.append(F['fx_p'])
        fy.append(F['fy_p'])
        fz.append(F['fz_p'])
        cx.append(F['cx'])
        cy.append(F['cy'])
        M_O.append(F['M_O'])
          
    for th in np.linspace(scroll_geo.theta_d(ScrollComp.geo)+0.001,2*pi,100):
        tt.append(th+2*pi)
        F = scroll_geo.D1_forces(th,geo=ScrollComp.geo,poly=True)
        fxp.append(F['fxp_poly'])
        fyp.append(F['fyp_poly'])
        fx.append(F['fx_p'])
        fy.append(F['fy_p'])
        fz.append(F['fz_p'])
        cx.append(F['cx'])
        cy.append(F['cy'])
        M_O.append(F['M_O'])
    
    import pylab
    pylab.plot(fx,fy,'-',cxp,fyp,'s',mfc='none')
    pylab.show()
    return
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
    Tc =  13 + 273.15
    Tin = Tc + 20
    DT_sc = 7
    pe = CP.Props('P','T',Te,'Q',1.0,Ref)
    pc = CP.Props('P','T',Tc,'Q',1.0,Ref)
    inletState = State.State(Ref,{'T':Tin,'P':pe})
    T2s = ScrollComp.isentropic_outlet_temp(inletState,pc)
    outletState = State.State(Ref,{'T':T2s,'P':pc})
    print 'pc',pc
    
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
                                        initialState = injState1,
                                        becomes = 'injCV.1'
                                        )
                          )
    
    ScrollComp.auto_add_CVs(inletState, outletState)
    ScrollComp.auto_add_leakage(flankFunc = ScrollComp.FlankLeakage, 
                                radialFunc = ScrollComp.RadialLeakage)
    
    ScrollComp.add_flow(FlowPath(key1='inlet.2', 
                                 key2='sa', 
                                 MdotFcn=ScrollComp.Inlet_sa,
                                 MdotFcn_kwargs = dict(X_d = 0.3)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1='sa', 
                                 key2='s1',
                                 MdotFcn=ScrollComp.SA_S,
                                 MdotFcn_kwargs = dict(X_d = 0.3)
                                 )
                        )
    ScrollComp.add_flow(FlowPath(key1='sa', 
                                 key2='s2', 
                                 MdotFcn=ScrollComp.SA_S,
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
        debug_plots(ScrollComp)

    print 'time taken',clock()-t1
    
    if Injection:
        print ScrollComp.FlowsProcessed.mean_mdot
        ha = CP.Props('H','T',Tc-DT_sc,'P',pc,'R404A')
        hb = CP.Props('H','T',Tsat,'Q',0.0,'R404A')
        hc = CP.Props('H','T',Tsat,'Q',1.0,'R404A')
        
        ScrollComp.injection_massflow_ratio = (ha-hb)/(hc-ha)
        print 'enthalpies',ha,hb,hc,'x',ScrollComp.injection_massflow_ratio
    
    debug_plots(ScrollComp)
    
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
            Compressor(0.3)
            FP = open('results.csv','w')
            for f in np.linspace(0.2,0.8,5):
                S = Compressor(f)
                x_sim = str(S.FlowsProcessed.mean_mdot['injection.1']/S.FlowsProcessed.mean_mdot['inlet.1']) 
                x_cycle = str(S.injection_massflow_ratio)
                FP.write('f, '+str(f)+','+str(S.FlowsProcessed.mean_mdot['injection.1'])+','+str(S.FlowsProcessed.mean_mdot['inlet.1'])+','+x_sim+','+x_cycle+'\n')
            FP.close()       
        else:
            Compressor()

