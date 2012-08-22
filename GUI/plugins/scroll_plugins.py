import sys; 
sys.path.append('..')

import pdsim_plugins
import pdsim_panels
from PDSim.flow.flow import FlowPath
from PDSim.core.core import Tube
from PDSim.core.containers import ControlVolume
from math import pi
import wx

class ScrollInjectionPlugin(pdsim_plugins.PDSimPlugin):
    
    short_description = 'Refrigerant injection for scroll'
    
    def __init__(self):
        """
        A plugin that adds the injection ports for the scroll compressor
        """
        self._activated = False
        
    def activate(self, event = None):
        #: The inputs toolbook that contains all the input panels
        ITB = self.GUI.MTB.InputsTB
        
        if not self._activated:
            #Add the panel to the inputs panel
            self.injection_panel = pdsim_panels.InjectionInputsPanel(ITB)
            ITB.AddPage(self.injection_panel,"Injection")
            self._activated = True
        else:
            page_names = [ITB.GetPageText(I) for I in range(ITB.GetPageCount())]
            I = page_names.index("Injection")
            ITB.RemovePage(I)
            self.injection_panel.Destroy()
            self._activated = False
            
    def apply(self, ScrollComp, **kwargs):
        """
        Add the necessary things for the scroll compressor injection
        
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
                                 L=L,
                                 ID=ID,
                                 mdot=0.001, 
                                 State1=ScrollComp.CVs['injCV.1'].State.copy(),
                                 fixed=1,
                                 TubeFcn=ScrollComp.TubeCode
                                 )
                            )
        ScrollComp.add_flow(FlowPath(key1='injection.2',
                                     key2='injCV.1',
                                     MdotFcn=ScrollComp.InjectionTubeFM))
        