
# If being run from the folder that contains the PDSim source tree, 
# remove the current location from the python path and use the 
# site-packages version of PDSim
import sys, os
current_path = os.path.abspath(os.curdir)
if current_path in sys.path and os.path.exists(os.path.join('PDSim','__init__.py')):
    i = sys.path.index(current_path)
    sys.path.pop(i)
else:
    print current_path,sys.path
    
from PDSim.flow.flow import FlowPath
from PDSim.scroll.core import Scroll
from PDSim.core.containers import ControlVolume
from PDSim.core.core import Tube
from PDSim.plot.plots import debug_plots
from CoolProp import State
from CoolProp import CoolProp as CP
from math import pi
    

        
def Compressor():
    
    ScrollComp=Scroll()
    #This runs if the module code is run directly
    ScrollComp.set_scroll_geo(104.8e-6,1.61,0.004,0.005) #Set the scroll wrap geometry
    ScrollComp.set_disc_geo('2Arc',r2='PMP')
    ScrollComp.geo.delta_flank=12e-6
    ScrollComp.geo.delta_radial=12e-6
    ScrollComp.omega=3500/60*2*pi
    ScrollComp.Tamb = 298.0
        
    Ref='Air'
#    State.debug(0)
#    CP.set_1phase_LUT_params(Ref,30,30,250,1000,180,5000)
#    CP.UseSinglePhaseLUT(True)
#    State.set_1phase_LUT_params(Ref,30,30,250,1000,180,5000)
#    State.LUT(True)
    
    inletState = State.State(Ref,{'T':300,'P':310})
    outletState = State.State(Ref,{'T':400,'P':1200})
    
    mdot_guess = inletState.rho*ScrollComp.Vdisp*ScrollComp.omega/(2*pi)
    
    ScrollComp.add_tube(Tube(key1='inlet.1',key2='inlet.2',L=0.03,ID=0.02,
                             mdot=mdot_guess, State1=inletState.copy(),
                             fixed=1,TubeFcn=ScrollComp.TubeCode))
    ScrollComp.add_tube(Tube(key1='outlet.1',key2='outlet.2',L=0.03,ID=0.02,
                             mdot=mdot_guess, State2=outletState.copy(),
                             fixed=2,TubeFcn=ScrollComp.TubeCode))
    
    ScrollComp.auto_add_CVs(inletState, outletState)
    ScrollComp.auto_add_leakage(flankFunc = ScrollComp.FlankLeakage, 
                                radialFunc = ScrollComp.RadialLeakage)
    
    print "still missing radial leakage terms"
    
    ScrollComp.add_flow(FlowPath(key1='inlet.2',key2='sa',MdotFcn=ScrollComp.Inlet_sa))
    ScrollComp.add_flow(FlowPath(key1='sa',key2='s1',MdotFcn=ScrollComp.SA_S))
    ScrollComp.add_flow(FlowPath(key1='sa',key2='s2',MdotFcn=ScrollComp.SA_S))
    
    ScrollComp.add_flow(FlowPath(key1='outlet.1',key2='ddd',MdotFcn=ScrollComp.Discharge))
    ScrollComp.add_flow(FlowPath(key1='outlet.1',key2='dd',MdotFcn=ScrollComp.Discharge))
    ScrollComp.add_flow(FlowPath(key1='d1',key2='dd',MdotFcn=ScrollComp.D_to_DD))
    ScrollComp.add_flow(FlowPath(key1='d2',key2='dd',MdotFcn=ScrollComp.D_to_DD))

    from time import clock
    t1=clock()
    ScrollComp.RK45_eps = 1e-6
    ScrollComp.EulerN = 10000
    ScrollComp.solve(key_inlet='inlet.1',key_outlet='outlet.2',
                 step_callback=ScrollComp.step_callback, 
                 endcycle_callback=ScrollComp.endcycle_callback,
                 heat_transfer_callback=ScrollComp.heat_transfer_callback,
                 lump_energy_balance_callback=ScrollComp.lump_energy_balance_callback,
                 solver_method='RK45',
                 hmin=2*pi/(100000),
                 UseNR = False,
                 OneCycle = True
                 )
    print 'time taken',clock()-t1
    
    debug_plots(ScrollComp)

    
if __name__=='__main__':
        
    profile=False
    if profile==True:
        import line_profiler as LP
        profiler=LP.LineProfiler(FlowPath.calculate,Scroll.CVs.updateStates)
        profiler.run("Compressor()")
        profiler.print_stats()
    else:
        Compressor()

