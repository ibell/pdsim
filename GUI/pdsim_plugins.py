from PDSim.flow.flow import FlowPath
from PDSim.core.core import Tube
from PDSim.core.containers import ControlVolume
from math import pi

class PDSimPlugin(object):
    def __init__(self):
        """
        A base class that is to be subclassed for other plugins to PDSim.
        
        Plugins are structured in order to be as flexible as possible.
        """
        pass
    
    def apply(self, sim, **kwargs):
        """
        Apply the plugin's code - it can do whatever it wants to the simulation 
        """
        raise NotImplementedError("Subclasses of PDSimPlugin must provide the apply() function")
    
    def validate(self):
        pass
    
    
class ScrollInjectionPlugin(PDSimPlugin):
    
    def __init__(self, injection_panel):
        """
        A plugin that adds the injection ports for the scroll compressor
        
        Parameters
        ----------
        injection_panel
        
        """
        self.injection_panel = injection_panel
    
        
    def apply(self, ScrollComp, **kwargs):
        """
        Add the necessary things
        
        Parameters
        ----------
        ScrollComp : Scroll instance
        """
        phi = float(self.injection_panel.InjectionElement.phival.GetValue())
        L = float(self.injection_panel.InjectionElement.Lval.GetValue())
        ID = float(self.injection_panel.InjectionElement.IDval.GetValue())
        injState1 = self.injection_panel.InjectionElement.state.GetState().copy()
        V_tube = L*pi*ID**2/4.0
        ScrollComp.add_CV(ControlVolume(key ='injCV.1',
                                        VdVFcn = ScrollComp.V_injection,
                                        VdVFcn_kwargs = dict(V_tube = V_tube),
                                        initialState = injState1,
                                        becomes = 'injCV.1'
                                        )
                          )
        
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
                                 mdot=0.001, 
                                 State1=ScrollComp.CVs['injCV.1'].State.copy(),
                                 fixed=1,
                                 TubeFcn=ScrollComp.TubeCode
                                 )
                            )
        ScrollComp.add_flow(FlowPath(key1='injection.2',
                                     key2='injCV.1',
                                     MdotFcn=ScrollComp.InjectionTubeFM))