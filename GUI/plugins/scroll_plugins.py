import sys
sys.path.append('..')

import pdsim_plugins
import pdsim_panels
from PDSim.flow.flow import FlowPath
from PDSim.core.core import Tube
from PDSim.core.containers import ControlVolume
from PDSim.scroll import scroll_geo
from math import pi
import wx
from PDSim.scroll.plots import plotScrollSet

LabeledItem = pdsim_panels.LabeledItem

class InjectionViewerDialog(wx.Dialog):
    """
    A dialog with simple plot of the locations of the injection ports
    """
    def __init__(self, theta, geo):
        wx.Dialog.__init__(self, parent = None)
        
        self.theta = theta
        
        self.geo = geo
        self.PP = pdsim_panels.PlotPanel(self)
        self.ax = self.PP.figure.add_axes((0,0,1,1))
        
        #Plot the injection ports (symmetric)
        plotScrollSet(theta, geo, axis = self.ax)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.PP)
        sizer.Layout()
        self.Fit()
        
    def add_port(self, phi, inner_outer):
        scroll_geo.plot_injection_ports(self.theta, self.geo, phi, self.ax, inner_outer)
        
class InjectionPortPanel(wx.Panel):
    """
    A panel with the values for one injection port
    """
    def __init__(self, parent):
        wx.Panel.__init__(self,parent)
        
        self.AddPort = wx.Button(self, label='+',style = wx.ID_REMOVE)
        self.RemovePort = wx.Button(self,label='-',style = wx.ID_REMOVE)
        self.AddPort.Bind(wx.EVT_BUTTON,lambda(event):self.Parent.OnAddPort(self))
        self.RemovePort.Bind(wx.EVT_BUTTON,lambda(event):self.Parent.OnRemovePort(self))
        
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        button_sizer.Add(self.AddPort,1,wx.EXPAND)
        button_sizer.Add(self.RemovePort)
        
        element_sizer = wx.FlexGridSizer(cols=2)
        element_sizer.Add(wx.StaticText(self,label="Involute Angle"))
        self.phi_inj_port = wx.TextCtrl(self,value="7.141") 
        self.phi_inj_port.SetToolTipString('If you want symmetric injection ports, the involute angle on the inner involute should be pi radians greater than that on the outer involute.  View the ports to be sure')
        element_sizer.Add(self.phi_inj_port)
        element_sizer.Add(wx.StaticText(self,label="Neighbor Involute"))
        self.involute = wx.ComboBox(self)
        self.involute.AppendItems(['Outer involute','Inner involute'])
        self.involute.SetSelection(0)
        self.involute.SetEditable(False)
        element_sizer.Add(self.involute)
        self.check_valve = wx.CheckBox(self,label="Uses check valve")
        element_sizer.Add(self.check_valve)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(button_sizer)
        sizer.Add(element_sizer)
        self.SetSizer(sizer)
        sizer.Layout()
        
    def get_values(self):
        """
        Returns a tuple of involute_angle, inner_outer, check_valve
        """
        if self.involute.GetStringSelection() == 'Outer involute':
            inner_outer = 'o'
        elif self.involute.GetStringSelection() == 'Inner involute':
            inner_outer = 'i'
        else:
            raise ValueError
        
        return float(self.phi_inj_port.GetValue()), inner_outer, self.check_valve.IsChecked()
        
        
class InjectionElementPanel(wx.Panel):
    """
    A panel with the injection values for one injection line with a box around it
    """
    def __init__(self, parent,index):
        wx.Panel.__init__(self,parent)
        
        #Inputs Toolbook
        ITB = self.GetTopLevelParent().MTB.InputsTB
        Fluid = None
        for panel in ITB.panels:
            if panel.Name == 'StatePanel':
                Fluid = panel.SuctionState.GetState().Fluid
                break
        if Fluid is None:
            raise ValueError('StatePanel not found in Inputs Toolbook')
        
        #You can only inject the same refrigerant as at the suction so fix the fluid
        self.state = pdsim_panels.StatePanel(self, Fluid=Fluid, Fluid_fixed = True)
        
        self.Llabel,self.Lval = LabeledItem(self, label='Length of injection line',value='1.0')
        self.IDlabel,self.IDval = LabeledItem(self, label='Inner diameter of injection line',value='0.01')
        line_sizer = wx.FlexGridSizer(cols = 2)
        line_sizer.AddMany([self.Llabel,self.Lval])
        line_sizer.AddMany([self.IDlabel,self.IDval])
        
        self.RemoveButton = wx.Button(self,label='Remove Line')
        self.RemoveButton.Bind(wx.EVT_BUTTON, lambda event: self.Parent.RemoveInjection(self))
        
        self.SizerBox = wx.StaticBox(self, label = "Injection line #"+str(index))
        self.SBSSizer = wx.StaticBoxSizer(self.SizerBox, wx.VERTICAL)
        
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.RemoveButton)
        
        text = wx.StaticText(self, label='Parameters:')
        font = text.GetFont()
        font.SetUnderlined(True)
        text.SetFont(font)
        sizer.Add(text)
        sizer.AddSpacer(10)
        sizer.Add(line_sizer)
        sizer.AddSpacer(10)
        text = wx.StaticText(self, label='State:')
        font = text.GetFont()
        font.SetUnderlined(True)
        text.SetFont(font)
        sizer.Add(text)
        sizer.Add(self.state)
        text = wx.StaticText(self, label='Injection Ports:')
        font = text.GetFont()
        font.SetUnderlined(True)
        text.SetFont(font)
        sizer.Add(text)
        sizer.Add(InjectionPortPanel(self))
        self.InnerSizer = sizer
        self.SBSSizer.Add(sizer)
        self.SBSSizer.Layout()
        self.SetSizer(self.SBSSizer)
        self.Nports = 1
        
    def OnRemovePort(self,port):
        if self.Nports > 1:
            port.Destroy()
            self.SBSSizer.Layout()
            self.Refresh()
            self.Nports -= 1
            self.Fit()
        
    def OnAddPort(self, event):
        self.SBSSizer.Add(InjectionPortPanel(self))
        self.SBSSizer.Layout()
        self.Nports += 1
        self.GetSizer().Layout()
        self.Refresh()
        self.Fit()
        
class InjectionInputsPanel(pdsim_panels.PDPanel):
    """
    The container panel for all the injection ports and injection data 
    """ 
    def __init__(self, parent):
        pdsim_panels.PDPanel.__init__(self,parent)
        
        #Add the header row of buttons
        self.View = wx.Button(self, label='View')
        self.View.Bind(wx.EVT_BUTTON, self.OnView)
        self.AddInjection = wx.Button(self, label='Add Injection Line')
        self.AddInjection.Bind(wx.EVT_BUTTON, self.OnAddInjection)
        buttons_sizer = wx.BoxSizer(wx.HORIZONTAL)
        buttons_sizer.Add(self.AddInjection)
        buttons_sizer.Add(self.View)
        
        sizer = wx.FlexGridSizer(cols = 1)
        sizer.Add(buttons_sizer)
        sizer.AddSpacer(10)
        sizer.Layout()
        self.SetSizer(sizer)
        self.Nterms = 0
        
    def OnAddInjection(self, event):
        """
        Add an injection line to the injection panel
        """
        IE = InjectionElementPanel(self,self.Nterms+1)
        self.GetSizer().Add(IE)
        self.Nterms += 1
        self.GetSizer().Layout()
        self.Refresh()
        
    def RemoveInjection(self, injection):
        """
        Remove the given injection term
        """
        injection.Destroy()
        self.Nterms -= 1
        #Renumber the injection panels
        I=1
        for child in self.Children:
            if isinstance(child,InjectionElementPanel):
                child.SizerBox.SetLabel("Injection line #"+str(I))
                I+=1
        self.GetSizer().Layout()
        self.Refresh()
        
    def OnView(self, event):
        geo = self.GetTopLevelParent().MTB.InputsTB.panels[0].Scroll.geo
        dlg = InjectionViewerDialog(pi/2,geo)
        
        #IEPs are children that are instances of InjectionElementPanel class
        IEPs = [child for child in self.Children if isinstance(child,InjectionElementPanel)]
        for IEP in IEPs:
            for child in IEP.Children:
                if isinstance(child,InjectionPortPanel):
                    phi,inner_outer,check_valve = child.get_values()
                    dlg.add_port(phi, inner_outer)
                    
        dlg.ShowModal()
        dlg.Destroy()
        
        
        
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
            self.injection_panel = InjectionInputsPanel(ITB)
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
        
        #IEPs are children of injection_panel that are instances of InjectionElementPanel class
        IEPs = [child for child in self.injection_panel.Children if isinstance(child,InjectionElementPanel)]
        for i,IEP in enumerate(IEPs):
            L = float(IEP.Lval.GetValue())
            ID = float(IEP.IDval.GetValue())
            injState = IEP.state.GetState().copy()
            V_tube = L*pi*ID**2/4.0
            
            CVkey = 'injCV.'+str(i+1)
            #Add the control volume for the injection line
            ScrollComp.add_CV(ControlVolume(key = CVkey,
                                            VdVFcn = ScrollComp.V_injection,
                                            VdVFcn_kwargs = dict(V_tube = V_tube),
                                            initialState = injState,
                                            )
                              )
            
            #Add the tube for the injection line
            ScrollComp.add_tube(Tube(key1='injection_line.'+str(i+1)+'.1',
                                     key2='injection_line.'+str(i+1)+'.2',
                                     L=L,
                                     ID=ID,
                                     mdot=0.001, 
                                     State1=ScrollComp.CVs[CVkey].State.copy(),
                                     fixed=1,
                                     TubeFcn=ScrollComp.TubeCode
                                     )
                                )
            
            #Add the flow model between the injection line tube and the injection CV 
            ScrollComp.add_flow(FlowPath(key1='injection_line.'+str(i+1)+'.2',
                                         key2=CVkey,
                                         MdotFcn=ScrollComp.IsentropicNozzleFM,
                                         MdotFcn_kwargs = dict(A = pi*ID**2/4)
                                         )
                                )
            
            for child in IEP.Children:
                if isinstance(child,InjectionPortPanel):
                    phi,inner_outer,check_valve = child.get_values()
                    
                    #Figure out which CV are in contact with this location for the injection port
                    partner_key_start = ScrollComp._get_injection_CVkey(phi, 0*pi, inner_outer)
                    partner_key_end = ScrollComp._get_injection_CVkey(phi, 2*pi, inner_outer)
                    
                    #Add the CV that start and end the rotation connected to the port
                    for partner_key in [partner_key_start, partner_key_end]:
                        #Injection flow paths
                        ScrollComp.add_flow(FlowPath(key1= partner_key, 
                                                     key2 = CVkey, 
                                                     MdotFcn=ScrollComp.Injection_to_Comp,
                                                     MdotFcn_kwargs = dict(phi = phi,
                                                                           inner_outer = inner_outer,
                                                                           check_valve = check_valve)
                                                    )
                                            )
        
#    def output_dictionary(self):
#        """
#        Return a dictionary of column options for the plugin that should be
#        added to the output on the OutputDataPanel
#        
#        This function is called before the simulation is run 
#         
#        Returns
#        -------
#        output_dict : dictionary
#            A dict with keys that are attributes of the ScrollComp instance, and
#            values that are the human-readable string representation of this term
#        """
#        return {'':
#                }
        
        