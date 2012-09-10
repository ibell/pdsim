# -*- coding: utf-8 -*-
import pdsim_panels
import wx
import numpy as np
from math import pi
from PDSim.scroll.core import Scroll
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as WXCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as WXToolbar
from PDSim.scroll.plots import plotScrollSet, ScrollAnimForm
from PDSim.flow.flow import FlowPath
from PDSim.core.core import Tube
from CoolProp import State as CPState

LabeledItem = pdsim_panels.LabeledItem

class PlotPanel(wx.Panel):
    def __init__(self, parent, **kwargs):
        wx.Panel.__init__(self, parent, size = (300,200), **kwargs)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.figure = mpl.figure.Figure(dpi=100, figsize=(2, 2))
#        self.figure.set_figwidth(2.0)
#        self.figure.set_figheight(2.0)
        self.canvas = WXCanvas(self, -1, self.figure)
#        self.canvas.resize(200,200)
        #self.toolbar = WXToolbar(self.canvas)
        #self.toolbar.Realize()
        sizer.Add(self.canvas)
        #sizer.Add(self.toolbar)
        self.SetSizer(sizer)
        sizer.Layout()
        
from wx.lib.scrolledpanel import ScrolledPanel
class GeometryPanel(pdsim_panels.PDPanel):
    """
    The geometry panel of the scroll compressor
    Loads all parameters from the configuration file
    """
    def __init__(self,parent,configfile,**kwargs):
        pdsim_panels.PDPanel.__init__(self,parent,**kwargs)
        
        #Now we are going to put everything into a scrolled window
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        scrolled_panel = ScrolledPanel(self, size=(-1,-1),
                                 style = wx.TAB_TRAVERSAL, name="panel1")
        scrolled_panel.SetScrollbars(1,1,1,1)
        
        #Loads all the parameters from the config file
        configdict, descdict = self.get_from_configfile('GeometryPanel')
        
        # Things in self.items are linked through to the module code where 
        # it attempts to set the attribute.  They are also automatically
        # written to configuration file
        self.items = [
        dict(attr='Vdisp'),
        dict(attr='Vratio'),
        dict(attr='t'),
        dict(attr='ro'),
        dict(attr='phi_fi0'),
        dict(attr='phi_fis'),
        dict(attr='phi_fos'),
        dict(attr='delta_flank'),
        dict(attr='delta_radial'),
        ]
        
        sizerInputs = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        
        self.ConstructItems(self.items, sizerInputs, configdict, descdict, parent = scrolled_panel)
        
        self.UseOffsetLabel = wx.StaticText(scrolled_panel, label = 'Use Offset scrolls')
        self.UseOffset = wx.CheckBox(scrolled_panel)
        self.UseOffset.SetValue(self._use_offset)
        del self._use_offset
        self.UseOffset.Bind(wx.EVT_CHECKBOX, self.OnChangeOffset)
        sizerInputs.AddMany([self.UseOffsetLabel,self.UseOffset])
        
        items2 = [dict(attr = 'delta_offset')]
        self.ConstructItems(items2, sizerInputs, configdict, descdict, parent = scrolled_panel)
        self.items += items2
        
        for item in self.items:
            setattr(self,item['attr'],item['textbox'])
        
        kwargs = dict(label = u"\u03D5_i0 [radian]",
                      tooltip = 'Initial involute angle for inner involute'
                      )
        self.phi_i0_label, self.phi_i0 = LabeledItem(scrolled_panel, **kwargs)
        
        self.phi_is_label, self.phi_is= LabeledItem(scrolled_panel,
                                                       label=u"\u03D5_is [radian]")

        width = max([item['label'].GetEffectiveMinSize()[0] for item in self.items])
        self.phi_is_label.SetMinSize((width,-1))
        self.phi_ie_label, self.phi_ie= LabeledItem(scrolled_panel,
                                                       label=u"\u03D5_ie [radian]")
        self.phi_o0_label, self.phi_o0= LabeledItem(scrolled_panel,
                                                       label=u"\u03D5_o0 [radian]")
        self.phi_os_label, self.phi_os= LabeledItem(scrolled_panel,
                                                       label=u"\u03D5_os [radian]")
        self.phi_oe_label, self.phi_oe= LabeledItem(scrolled_panel,
                                                       label=u"\u03D5_oe [radian]")
        self.rb_label, self.rb = LabeledItem(scrolled_panel,
                                             label="rb [m]")
        self.hs_label, self.hs= LabeledItem(scrolled_panel,
                                            label="hs [m]")
        
        self.phi_i0.Enable(False)
        self.phi_is.Enable(False)
        self.phi_ie.Enable(False)
        self.phi_o0.Enable(False)
        self.phi_os.Enable(False)
        self.phi_oe.Enable(False)
        self.rb.Enable(False)
        self.hs.Enable(False)
        
        self.PP = PlotPanel(scrolled_panel)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(self.PP,0,wx.EXPAND)
        anibutton = wx.Button(scrolled_panel, label = 'Animate')
        anibutton.Bind(wx.EVT_BUTTON,self.OnAnimate)
        hsizer.Add(anibutton,1,wx.EXPAND)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizerInputs)
        sizer.AddSpacer(10)
        sizer.Add(hsizer)
        sizer.AddSpacer(10)
        
        self.ax = self.PP.figure.add_axes((0,0,1,1))

        fgsGeoAnglesInputs = wx.FlexGridSizer(cols = 2)
        fgsGeoAnglesInputs.AddMany([
                     self.phi_i0_label,self.phi_i0, 
                     self.phi_is_label,self.phi_is,
                     self.phi_ie_label,self.phi_ie,
                     self.phi_o0_label,self.phi_o0, 
                     self.phi_os_label,self.phi_os,
                     self.phi_oe_label,self.phi_oe,
                     self.rb_label,self.rb,
                     self.hs_label,self.hs
                     ])
        sizer.Add(fgsGeoAnglesInputs)

        
        for item in self.items:
            item['textbox'].Bind(wx.EVT_KILL_FOCUS, self.OnChangeParam)
        self.UseOffset.Bind(wx.EVT_CHECKBOX,self.OnChangeParam)
        delta_offset_item = self._get_item_by_attr('delta_offset')['textbox']
        delta_offset_item.Bind(wx.EVT_KILL_FOCUS, self.OnChangeParam)
        
        self.items += items2
        # Keep a local copy of the scroll in order to be able to use the 
        # set_scroll_geo and set_disc_geo functions
        self.Scroll=Scroll()
        
        #Do the layout of all the panels
        scrolled_panel.SetSizer(sizer)
        main_sizer.Add(scrolled_panel,1,wx.EXPAND)
        self.SetSizer(main_sizer)   
        
        self.OnChangeOffset()
        self.OnChangeParam()
          
    def post_prep_for_configfile(self):
        return 'use_offset = bool, ' + str(self.UseOffset.IsChecked()) + '\n'
    
    def post_get_from_configfile(self,k,v):
        
        if k == 'use_offset':
            self._use_offset = (v.split(',')[1].strip().lower() == 'true')
    
    def skip_list(self):
        """
        Returns a list of atttributes to skip setting in set_params() function
        from the base class PDPanel
        """
        return ['Vdisp','Vratio','t','ro']
        
    def OnChangeOffset(self, event = None):
        delta_offset_item = self._get_item_by_attr('delta_offset')['textbox']
        if self.UseOffset.IsChecked():
            delta_offset_item.Enable(True)
        else:
            delta_offset_item.Enable(False)
        
    def OnAnimate(self, event = None):
        SAF = ScrollAnimForm(self.Scroll.geo)
        SAF.Show()
        
    def OnChangeParam(self, event = None):
        if event is not None:
            event.Skip()
        Vdisp=float(self.Vdisp.GetValue())
        Vratio=float(self.Vratio.GetValue())
        t=float(self.t.GetValue())
        ro=float(self.ro.GetValue())
        phi_i0 = float(self.phi_fi0.GetValue())
        phi_is = float(self.phi_fis.GetValue())
        phi_os = float(self.phi_fos.GetValue())
        
        #Set the scroll wrap geometry
        self.Scroll.set_scroll_geo(Vdisp,Vratio,t,ro,
                                   phi_i0 = phi_i0, 
                                   phi_os = phi_os,
                                   phi_is = phi_is,) 
        self.Scroll.set_disc_geo('2Arc', r2='PMP')
        
        if self.UseOffset.IsChecked():
            self.Scroll.geo.phi_ie_offset = pi
            delta_offset_item = self._get_item_by_attr('delta_offset')['textbox']
            self.Scroll.geo.delta_suction_offset = float(delta_offset_item.GetValue())
        else:
            self.Scroll.geo.phi_ie_offset = 0
            self.Scroll.geo.delta_suction_offset = 0.0
        
        self.phi_i0.SetValue(str(self.Scroll.geo.phi_i0))
        self.phi_is.SetValue(str(self.Scroll.geo.phi_is))
        self.phi_ie.SetValue(str(self.Scroll.geo.phi_ie))
        self.phi_o0.SetValue(str(self.Scroll.geo.phi_o0))
        self.phi_os.SetValue(str(self.Scroll.geo.phi_os))
        self.phi_oe.SetValue(str(self.Scroll.geo.phi_oe))
        self.rb.SetValue(str(self.Scroll.geo.rb))
        self.hs.SetValue(str(self.Scroll.geo.h))
        
        self.ax.cla()

        plotScrollSet(pi/4.0, 
                      axis = self.ax, 
                      geo = self.Scroll.geo,
                      offsetScroll = self.Scroll.geo.phi_ie_offset > 0)
        self.PP.canvas.draw()
    
    def post_set_params(self, scroll):
        Vdisp=float(self.Vdisp.GetValue())
        Vratio=float(self.Vratio.GetValue())
        t=float(self.t.GetValue())
        ro=float(self.ro.GetValue())
        
        scroll.set_scroll_geo(Vdisp,Vratio,t,ro) #Set the scroll wrap geometry
        scroll.set_disc_geo('2Arc', r2='PMP')
        scroll.geo.delta_flank = float(self.delta_flank.GetValue())
        scroll.geo.delta_radial = float(self.delta_radial.GetValue())
        
        if self.UseOffset.IsChecked():
            scroll.geo.phi_ie_offset = pi
        else:
            scroll.geo.phi_ie_offset = 0
    
    
    def collect_output_terms(self):
        """
        
        Hacked to change the parameters in the output table
        """
        print 'changing params'
        Main = self.GetTopLevelParent()
        Main.MTB.OutputsTB.DataPanel.change_output_attrs(dict(t = 'geo.t',
                                                              ro = 'geo.ro',
                                                              )
                                                         )
        return []
        
class FlankLeakageFlowChoice(pdsim_panels.MassFlowOption):
    
    def __init__(self,parent,**kwargs):
        pdsim_panels.MassFlowOption.__init__(self,parent,**kwargs)
    
    def model_options(self):
        return [
                dict(desc = 'Hybrid leakage model',
                     function_name = 'FlankLeakage',)
                ]
        
class RadialLeakageFlowChoice(pdsim_panels.MassFlowOption):
    
    def __init__(self,parent,**kwargs):
        pdsim_panels.MassFlowOption.__init__(self,parent,**kwargs)
    
    def model_options(self):
        return [
                dict(desc = 'Hybrid leakage model',function_name = 'RadialLeakage')
                ]

class SuctionFlowChoice(pdsim_panels.MassFlowOption):
    
    def __init__(self,parent,**kwargs):
        pdsim_panels.MassFlowOption.__init__(self,parent,**kwargs)
    
    def model_options(self):
        return [
                dict(desc = 'Isentropic nozzle',
                     function_name = 'SA_S',
                     params = [dict(attr = 'X_d',
                                    value = 1.0,
                                    desc = 'Tuning factor')]
                     )
                ]
                
class InletFlowChoice(pdsim_panels.MassFlowOption):
    
    def __init__(self,parent,**kwargs):
        pdsim_panels.MassFlowOption.__init__(self,parent,**kwargs)
        self.parent = parent
        
    def model_options(self):
        return [dict(desc = 'Isentropic nozzle',
                     function_name = 'IsentropicNozzleFM',
                     params = [dict(attr = 'X_d',
                                    value = 1.0,
                                    desc = 'Tuning factor')
                               ]
                     )
                ]
                
class MassFlowPanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, configfile,**kwargs):
    
        pdsim_panels.PDPanel.__init__(self, parent,**kwargs)
        
        #Loads all the parameters from the config file
        self.configdict,self.descdict = self.get_from_configfile('MassFlowPanel')
        
        self.items1 = [
        dict(attr='d_discharge'),
        dict(attr='inlet_tube_length'),
        dict(attr='inlet_tube_ID'),
        dict(attr='outlet_tube_length'),
        dict(attr='outlet_tube_ID')
        ]
        box_sizer = wx.BoxSizer(wx.VERTICAL)
        box_sizer.Add(wx.StaticText(self,-1,"Required Inputs"))
        box_sizer.Add(wx.StaticLine(self,-1,(25, 50), (300,1)))
        
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        self.ConstructItems(self.items1,sizer,self.configdict,self.descdict)

        box_sizer.Add(sizer)  
        
        box_sizer.AddSpacer(10)
        box_sizer.Add(wx.StaticText(self,-1,"Flow Models"))
        box_sizer.Add(wx.StaticLine(self,-1,(25, 50), (300,1)))
#        self.flankflow = FlankLeakageFlowChoice(parent = self,
#                                           label = 'Flank leakage')
#        self.radialflow = RadialLeakageFlowChoice(parent = self,
#                                           label = 'Radial leakage')
        self.suctionflow1 = SuctionFlowChoice(parent = self,
                                              key1 = 'sa',
                                              key2 = 's1',
                                              label = 'Flow to suction chamber #1')
        self.suctionflow2 = SuctionFlowChoice(parent = self,
                                              key1 = 'sa',
                                              key2 = 's2',
                                              label = 'Flow to suction chamber #2')    
        self.inletflow = InletFlowChoice(parent = self,
                                         key1 = 'inlet.2',
                                         key2 = 'sa',
                                         label = 'Flow into shell',
                                         )
        
        self.flows = [self.suctionflow1, self.suctionflow2, self.inletflow]
#        self.flows = [self.flankflow,self.radialflow,self.suctionflow,self.inletflow]
        
        box_sizer.AddMany(self.flows)
        
        self.resize_flows(self.flows)
        
        self.SetSizer(box_sizer)
        box_sizer.Layout()
        
        self.items=self.items1
        
    def resize_flows(self, flows):
        min_width = max([flow.label.GetSize()[0] for flow in flows])
        for flow in flows:
            flow.label.SetMinSize((min_width,-1))
            
    def post_set_params(self, simulation):
        #Create and add each of the flow paths based on the flow model selection
        for flow in self.flows:
            func_name = flow.get_function_name()
            #func is a pointer to the actual function in the simulation instance
            func = getattr(simulation, func_name)
            
            param_dict = {p['attr']:p['value'] for p in flow.params_dict}
            if isinstance(flow, InletFlowChoice):
                #Get the diameter from the inlet_tube_ID item
                D = float(self._get_item_by_attr('inlet_tube_ID')['textbox'].GetValue())
                #Apply the correction factor
                A = pi*D**2/4 * param_dict.pop('X_d')
                param_dict['A']=A
                
            simulation.add_flow(FlowPath(key1 = flow.key1,
                                         key2 = flow.key2,
                                         MdotFcn = func,
                                         MdotFcn_kwargs = param_dict
                                         )
                                )
            
        D = float(self._get_item_by_attr('d_discharge')['textbox'].GetValue())
        A_discharge_port = pi*D**2/4
        
        for _key in ['dd','ddd']:
            simulation.add_flow(FlowPath(key1='outlet.1', 
                                         key2=_key, 
                                         MdotFcn=simulation.IsentropicNozzleFM,
                                         MdotFcn_kwargs = dict(A = A_discharge_port)
                                         )
                                )
    
        if callable(simulation.Vdisp):
            Vdisp = simulation.Vdisp()
        else:
            Vdisp = simulation.Vdisp
        
        # Set omega and inlet state
        # omega is set in the conventional way, so you want to make sure you 
        # don't overwrite it with the GUI value if you are doing a parametric table
        # state is set using the special method for additional parametric terms
        
        parent = self.GetParent() #InputsToolBook
        for child in parent.GetChildren():
            if hasattr(child,'Name') and child.Name == 'StatePanel':
                if hasattr(simulation,'omega'): 
                    omega = simulation.omega
                    child.set_params(simulation)
                    child.post_set_params(simulation)
                    simulation.omega = omega
                else:
                    child.set_params(simulation)
                    child.post_set_params(simulation)
                    omega = simulation.omega
                        
        Vdot = Vdisp*simulation.omega/(2*pi)
        
        T2s = simulation.guess_outlet_temp(simulation.inletState, simulation.discharge_pressure)
        outletState=CPState.State(simulation.inletState.Fluid,{'T':T2s,'P':simulation.discharge_pressure})
        
        simulation.auto_add_leakage(flankFunc = simulation.FlankLeakage, 
                                    radialFunc = simulation.RadialLeakage)
            
        #Create and add each of the inlet and outlet tubes
        simulation.add_tube( Tube(key1='inlet.1',
                                  key2='inlet.2',
                                  L=simulation.inlet_tube_length, 
                                  ID=simulation.inlet_tube_ID,
                                  mdot=simulation.inletState.copy().rho*Vdot, 
                                  State1=simulation.inletState.copy(),
                                  fixed=1, 
                                  TubeFcn=simulation.TubeCode) )
    
        simulation.add_tube( Tube(key1='outlet.1',
                                  key2='outlet.2',
                                  L=simulation.outlet_tube_length,
                                  ID=simulation.outlet_tube_ID,
                                  mdot=outletState.copy().rho*Vdot, 
                                  State2=outletState.copy(),
                                  fixed=2,
                                  TubeFcn=simulation.TubeCode) )
        
        
class MechanicalLossesPanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, configfile,**kwargs):
    
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        #Loads all the parameters from the config file (case-sensitive)
        self.configdict, self.descdict = self.get_from_configfile('MechanicalLossesPanel')
        
        self.items = [
        dict(attr='eta_motor'),
        dict(attr='h_shell'),
        dict(attr='A_shell'),
        dict(attr='Tamb'),
        dict(attr='mu_oil'),
        dict(attr='D_upper_bearing'),
        dict(attr='L_upper_bearing'),
        dict(attr='c_upper_bearing'),
        dict(attr='D_crank_bearing'),
        dict(attr='L_crank_bearing'),
        dict(attr='c_crank_bearing'),
        dict(attr='D_lower_bearing'),
        dict(attr='L_lower_bearing'),
        dict(attr='c_lower_bearing'),
        
        dict(attr='thrust_friction_coefficient')
        ]
        
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        
        self.ConstructItems(self.items,sizer,self.configdict,self.descdict)

        self.SetSizer(sizer)
        sizer.Layout()