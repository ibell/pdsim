# -*- coding: utf-8 -*-
import pdsim_panels
import wx
from wx.lib.mixins.listctrl import TextEditMixin
from wx.lib.scrolledpanel import ScrolledPanel
import numpy as np
from math import pi
from PDSim.scroll.core import Scroll
import matplotlib as mpl
import textwrap

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as WXCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as WXToolbar
from PDSim.scroll.plots import plotScrollSet, ScrollAnimForm
from PDSim.flow.flow import FlowPath
from PDSim.core.core import Tube
from PDSim.core.motor import Motor
from PDSim.misc.datatypes import AnnotatedValue
from CoolProp import State as CPState
from pdsim_panels import LaTeXImageMaker, MotorChoices, PlotPanel
from datatypes import HeaderStaticText

LabeledItem = pdsim_panels.LabeledItem
        
class ReadOnlyLaTeXLabel(wx.Panel):
    """
    A stub panel to allow for a LaTeX image with an additional caption for units
    """
    def __init__(self, LaTeX, parent, remaining_label = ''):
        wx.Panel.__init__(self, parent = parent)
        
        # Sizer
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        # The objects
        img = LaTeXImageMaker(LaTeX, parent = self)
        lab = wx.StaticText(self, label = remaining_label)
        
        # Layout
        sizer.AddMany([img,lab])
        self.SetSizer(sizer)
        sizer.Layout()
        
    def GetValue(self):
        return self.textbox.GetValue()
    
    def SetValue(self, value):
        self.textbox.SetValue(value)
        
class ScrollWrapAnglesFrame(wx.Frame):
    def __init__(self, geo):
        wx.Frame.__init__(self, None)
        
        panel = wx.Panel(self)
        
        # The sizer for all the outputs
        sizer_for_outputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        label1 = ReadOnlyLaTeXLabel('$\phi_{i0}$', parent = panel, remaining_label='[rad]')
        self.phi_i0 = wx.TextCtrl(panel)
        self.phi_i0.SetEditable(False)
        label2 = ReadOnlyLaTeXLabel('$\phi_{is}$', parent = panel, remaining_label='[rad]')
        self.phi_is = wx.TextCtrl(panel)
        self.phi_is.SetEditable(False)
        label3 = ReadOnlyLaTeXLabel('$\phi_{ie}$', parent = panel, remaining_label='[rad]')
        self.phi_ie = wx.TextCtrl(panel)
        self.phi_ie.SetEditable(False)
        label4 = ReadOnlyLaTeXLabel('$\phi_{o0}$', parent = panel, remaining_label='[rad]')
        self.phi_o0 = wx.TextCtrl(panel)
        self.phi_o0.SetEditable(False)
        label5 = ReadOnlyLaTeXLabel('$\phi_{os}$', parent = panel, remaining_label='[rad]')
        self.phi_os = wx.TextCtrl(panel)
        self.phi_os.SetEditable(False)
        label6 = ReadOnlyLaTeXLabel('$\phi_{oe}$', parent = panel, remaining_label='[rad]')
        self.phi_oe = wx.TextCtrl(panel)
        self.phi_oe.SetEditable(False)
        label7 = ReadOnlyLaTeXLabel('$r_b$', parent = panel, remaining_label='[m]')
        self.rb = wx.TextCtrl(panel)
        self.rb.SetEditable(False)
        label8 = ReadOnlyLaTeXLabel('$h_s$', parent = panel, remaining_label='[m]')
        self.hs = wx.TextCtrl(panel)
        self.hs.SetEditable(False)
        
        #Set the values of each of the boxes
        self.phi_i0.SetValue(str(geo.phi_i0))
        self.phi_is.SetValue(str(geo.phi_is))
        self.phi_ie.SetValue(str(geo.phi_ie))
        self.phi_o0.SetValue(str(geo.phi_o0))
        self.phi_os.SetValue(str(geo.phi_os))
        self.phi_oe.SetValue(str(geo.phi_oe))
        self.rb.SetValue(str(geo.rb))
        self.hs.SetValue(str(geo.h))
        
        # Add all the output objects to the sizer for the outputs
        sizer_for_outputs.AddMany([label1, self.phi_i0,
                                   label2, self.phi_is,
                                   label3, self.phi_ie,
                                   label4, self.phi_o0,
                                   label5, self.phi_os,
                                   label6, self.phi_oe,
                                   label7, self.rb,
                                   label8, self.hs])
        
        #Do the layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizer_for_outputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(wx.Button(panel, label='Close'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        panel.SetSizer(sizer)
        sizer.Layout()
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(panel)
        self.SetSizer(main_sizer)
        main_sizer.Layout()
        self.SetClientSize(main_sizer.GetMinSize())
    
class GeometryPanel(pdsim_panels.PDPanel):
    """
    The geometry panel of the scroll compressor
    Loads all parameters from the configuration file
    """
    
    # Maps from key in config file to description of the term  
    desc_map = dict(Vdisp = ('Displacement of the machine [m\xb3/rev]','m^3'),
                    Vratio = ('Built-in volume ratio [-]','-'),
                    t = ('Thickness of the scroll wrap [m]','m'),
                    ro = ('Orbiting radius [m]','m'),
                    phi_fi0 = ('Initial involute angle of the inner involute of the fixed scroll [rad]','rad'),
                    phi_fis = ('Starting involute angle of the inner involute of the fixed scroll [rad]','rad'),
                    phi_fos = ('Starting involute angle of the outer involute of the fixed scroll [rad]','rad'),
                    use_offset = ('Use offset geometry',''),
                    delta_offset = ('Offset gap width [m]','m'),
                    delta_flank = ('Flank gap width [m]','m'),
                    delta_radial = ('Radial gap width [m]' ,'m'),
                    d_discharge = ('Discharge port diameter [m]','m'),
                    inlet_tube_length = ('Inlet tube length [m]','m'),
                    inlet_tube_ID = ('Inlet tube inner diameter [m]','m'),
                    outlet_tube_length = ('Outlet tube length [m]','m'),
                    outlet_tube_ID = ('Outlet tube inner diameter [m]','m')
                    )
    
    def __init__(self, parent, config, **kwargs):
        """
        Parameters
        ----------
        parent : wx.Panel
            The parent of this panel
        config : dict
            The section of the configuration file pertaining to the geometry panel
        """
        # Instantiate the base class
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        # Now we are going to put everything into a scrolled window
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # The scrolled panel
        scrolled_panel = ScrolledPanel(self, size = (-1,-1), style = wx.TAB_TRAVERSAL, name="panel1")
        scrolled_panel.SetScrollbars(1, 1, 1, 1)
        
        # The list for all the annotated objects
        self.annotated_values = []
        
        # The sizer for all the objects
        sizer_for_wrap_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        annotated_values = []
        # Loop over the first group of inputs
        for key in ['Vdisp','Vratio','t','ro','phi_fi0','phi_fis','phi_fos',
                    'use_offset','delta_offset','delta_flank','delta_radial']:
            # Get the annotation and the units for the term 
            annotation, units = self.desc_map[key]
            # Add the annotated object to the list of objects
            annotated_values.append(AnnotatedValue(key, config[key], annotation, units))
            
        self.ScrollWrapAnglesButton = wx.Button(scrolled_panel,label='View Scroll Wrap Angles')
        self.ScrollWrapAnglesButton.Bind(wx.EVT_BUTTON,self.OnShowWrapGeo)
        
        # Build the items and return the list of annotated GUI objects
        annotated_GUI_objects = self.construct_items(annotated_values, 
                                                     sizer = sizer_for_wrap_inputs, 
                                                     parent = scrolled_panel)
        
        # The sizer for all the tube objects
        sizer_for_tube_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the tube inputs
        annotated_values = []
        for key in ['inlet_tube_length', 'inlet_tube_ID', 'outlet_tube_length', 'outlet_tube_ID']:
            # Get the annotation and the units for the term 
            annotation, units = self.desc_map[key]
            # Add the annotated object to the list of objects
            annotated_values.append(AnnotatedValue(key, config[key], annotation, units))
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                     sizer = sizer_for_tube_inputs,
                                                     parent = scrolled_panel)
        
        
        # ---------------------------------------------------------------------
        # Register terms in the GUI database
        self.main.register_GUI_objects(annotated_GUI_objects)
        
        # Link callback for refreshing of this panel with changing any input
        # parameter 
        for o in annotated_GUI_objects:
            o.GUI_location.Bind(wx.EVT_KILL_FOCUS, self.OnRefresh)
        
        # The plot of the scroll wraps
        self.PP = PlotPanel(scrolled_panel)
        self.ax = self.PP.figure.add_axes((0, 0, 1, 1))
        anibutton = wx.Button(scrolled_panel, label = 'Animate')
        anibutton.Bind(wx.EVT_BUTTON, self.OnAnimate)
        plotwrapssizer = wx.BoxSizer(wx.HORIZONTAL)
        plotwrapssizer.Add(self.PP, 1, wx.EXPAND)
        plotwrapssizer.Add(anibutton, 0, wx.EXPAND)

        # Layout the sizers
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(HeaderStaticText(scrolled_panel, 'Scroll Wrap Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(plotwrapssizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(self.ScrollWrapAnglesButton, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_wrap_inputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(HeaderStaticText(scrolled_panel, 'Discharge Port Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(HeaderStaticText(scrolled_panel, 'Tube Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_tube_inputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        
        # Do the layout of the scrolled panel
        scrolled_panel.SetSizer(sizer)
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND)
        self.SetSizer(main_sizer)

        # Create a scroll model instance to hold the geometry 
        self.Scroll = Scroll()
        
        # Refresh the panel
        self.OnRefresh()
        
    def OnShowWrapGeo(self, event = None):
        if event is not None: event.Skip()
        
        frm = ScrollWrapAnglesFrame(self.Scroll.geo)
        frm.Show()
        
    def OnAnimate(self, event = None):
        SAF = ScrollAnimForm(self.Scroll.geo, size=(400,400))
        SAF.Show()
        
    def OnRefresh(self, event = None):
        if event is not None: event.Skip()
        
        def get(key):
            # Compact code to get a parameter from the main database
            return self.main.get_GUI_object_value(key)
        
        # Set the scroll wrap geometry
        self.Scroll.set_scroll_geo(get('Vdisp'),
                                   get('Vratio'),
                                   get('t'),
                                   get('ro'),
                                   phi_i0 = get('phi_fi0'),
                                   phi_os = get('phi_fos'),
                                   phi_is = get('phi_fis')
                                   )
        self.Scroll.set_disc_geo('2Arc', r2 = 0)
        
        if get('use_offset'):
            self.Scroll.geo.phi_ie_offset = pi
            self.Scroll.geo.delta_suction_offset = get('delta_offset')
            self.main.get_GUI_object('delta_offset').GUI_location.Enable(True)
        else:
            self.Scroll.geo.phi_ie_offset = 0
            self.Scroll.geo.delta_suction_offset = 0.0
            self.main.get_GUI_object('delta_offset').GUI_location.Enable(False)
        
        
        
        self.ax.cla()

        plotScrollSet(pi/4.0, 
                      axis = self.ax, 
                      geo = self.Scroll.geo,
                      offsetScroll = self.Scroll.geo.phi_ie_offset > 0)
        
        # Plot the discharge port if the variable _d_discharge has been set
        try:
            d_discharge = get('d_discharge')
            t = np.linspace(0, 2*np.pi)
            x = self.Scroll.geo.xa_arc1 + d_discharge/2*np.cos(t)
            y = self.Scroll.geo.ya_arc1 + d_discharge/2*np.sin(t)
            self.ax.plot(x,y,'--')
        except KeyError:
            pass
            
        self.PP.canvas.draw()
        
    def get_script_chunks(self):
        
        def get(key):
            # Compact code to get a parameter from the main database
            return self.main.get_GUI_object(key).GetValue()
    
        if get('use_offset'):
            phi_ie_offset = str(pi)
        else:
            phi_ie_offset = str(0)
            
        #Parameters to be set in the string:
        str_params = dict(Vdisp = get('Vdisp'),
                          Vratio = get('Vratio'),
                          t = get('t'),
                          ro = get('ro'),
                          phi_i0 = get('phi_i0'),
                          phi_os = get('phi_os'),
                          phi_is = get('phi_is'),
                          phi_ie_offset = phi_ie_offset)

        return textwrap.dedent(
            """
            #Parameters from the GUI
            Vdisp = {Vdisp:s} #[m^3/rev]
            Vratio = {Vratio:s} #[-] 
            t = {t:s} #[m]
            ro = {ro:s} #[m]
            phi_i0 = {phi_i0:s} #[rad]
            phi_is = {phi_is:s} #[rad]
            phi_os = {phi_os:s} #[rad]
            
            #Set the scroll wrap geometry
            sim.set_scroll_geo(Vdisp, # Vdisp [m^3/rev]
                               Vratio, # Vratio [-]
                               t, # Thickness [m]
                               ro, # Orbiting radius [m]
                               phi_i0 = phi_i0, # [rad]
                               phi_os = phi_os, # [rad]
                               phi_is = phi_is) # [rad]
            sim.set_disc_geo('2Arc', r2 = 0) #hard-coded
            
            # self.delta_flank and self.delta_radial are set in the conventional way
            # using the parametric table
            sim.geo.delta_flank = sim.delta_flank
            sim.geo.delta_radial = sim.delta_radial
            
            sim.geo.phi_ie_offset = {phi_ie_offset:s}        
            """.format(**str_params))
    
    def collect_output_terms(self):
        """
        
        Hacked to change the parameters in the output table
        """
        print 'changing params'
        Main = self.GetTopLevelParent()
        Main.MTB.OutputsTB.DataPanel.change_output_attrs(dict(t = 'geo.t',
                                                              ro = 'geo.ro',))
        return []
        
class FlankLeakageFlowChoice(pdsim_panels.MassFlowOptionPanel):
    
    def __init__(self,parent,**kwargs):
        pdsim_panels.MassFlowOptionPanel.__init__(self,parent,**kwargs)
    
    def model_options(self):
        return [
                dict(desc = 'Hybrid leakage model',
                     function_name = 'FlankLeakage',)
                ]
        
class RadialLeakageFlowChoice(pdsim_panels.MassFlowOptionPanel):
    
    def __init__(self,parent,**kwargs):
        pdsim_panels.MassFlowOptionPanel.__init__(self,parent,**kwargs)
    
    def model_options(self):
        return [
                dict(desc = 'Hybrid leakage model',function_name = 'RadialLeakage')
                ]

class SuctionFlowChoice(pdsim_panels.MassFlowOptionPanel):
    
    def __init__(self,parent,**kwargs):
        pdsim_panels.MassFlowOptionPanel.__init__(self,parent,**kwargs)
    
    def model_options(self):
        return [
                dict(desc = 'Isentropic nozzle',
                     function_name = 'SA_S',
                     params = [dict(attr = 'X_d',
                                    value = 1.0,
                                    desc = 'Tuning factor')]
                     )
                ]
                
class InletFlowChoice(pdsim_panels.MassFlowOptionPanel):
    
    def __init__(self,parent,**kwargs):
        pdsim_panels.MassFlowOptionPanel.__init__(self,parent,**kwargs)
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

        self.suctionflow1 = SuctionFlowChoice(parent = self,
                                              key1 = 'sa',
                                              key2 = 's1',
                                              label = 'Flow to suction chamber #1')
        Xd_sa_s1,Xd_sa_s1_desc = self.get_from_configfile('MassFlowPanel','Xd_sa_s1')
        self.suctionflow1.set_attr('X_d', float(Xd_sa_s1))
        
        self.suctionflow2 = SuctionFlowChoice(parent = self,
                                              key1 = 'sa',
                                              key2 = 's2',
                                              label = 'Flow to suction chamber #2')    
        Xd_sa_s2,Xd_sa_s2_desc = self.get_from_configfile('MassFlowPanel','Xd_sa_s2')
        self.suctionflow2.set_attr('X_d', float(Xd_sa_s2))
        
        self.inletflow = InletFlowChoice(parent = self,
                                         key1 = 'inlet.2',
                                         key2 = 'sa',
                                         label = 'Flow into shell',
                                         )
        Xd_shell_sa,Xd_shell_sa_desc = self.get_from_configfile('MassFlowPanel','Xd_inlet.2_sa')
        self.inletflow.set_attr('X_d', float(Xd_shell_sa))
        
        self.flows = [self.suctionflow1, self.suctionflow2, self.inletflow]
        
        box_sizer.AddMany(self.flows)
        
        self.resize_flows(self.flows)
        
        self.SetSizer(box_sizer)
        box_sizer.Layout()
        
        self.items=self.items1
        
        # When discharge port diameter changes, update the value for the plot 
        # in the geometry panel
        self.d_discharge = self._get_item_by_attr('d_discharge')['textbox']
        self.d_discharge.Bind(wx.EVT_KILL_FOCUS, self.OnChangeDdisc) 
    
    def OnChangeDdisc(self, event = None):
        """
        Callback to set the variable _d_discharge in the Geometry Panel
        """
        GeoPanel = self.Parent.panels_dict['GeometryPanel']
        # Set the internal variable
        GeoPanel._d_discharge = float(self.d_discharge.GetValue())
        # Re-plot
        GeoPanel.OnChangeParam()
        
        
    def resize_flows(self, flows):
        """
        Resize the labels for the flows to all be the same size
        """
        min_width = max([flow.label.GetSize()[0] for flow in flows])
        for flow in flows:
            flow.label.SetMinSize((min_width,-1))
     
    def post_prep_for_configfile(self):
        """
        Prepare a string that saves the flow coefficients to file
        """
        s=''
        for flow in self.flows:
            param_dict = {p['attr']:p['value'] for p in flow.params_dict}
            s += 'Xd_{key1:s}_{key2:s} = float, {Xd:g}, flow coefficient\n'.format(key1 = flow.key1,
                                                          key2 = flow.key2,
                                                          Xd = param_dict.pop('X_d'))
        return s
        
    def get_script_chunks(self):
        
        Xd_dict = {}
        for flow in self.flows:
            param_dict = {p['attr']:p['value'] for p in flow.params_dict}
            Xd_dict[flow.key1+'-'+flow.key2] = str(param_dict.pop('X_d'))
          
          
        return textwrap.dedent(
            """
            # Add both the inlet and outlet tubes
            mdot_guess = inletState.rho*sim.Vdisp*sim.omega/(2*pi)
    
            sim.add_tube(Tube(key1 = 'inlet.1',
                              key2 = 'inlet.2',
                              L = sim.inlet_tube_length,
                              ID = sim.inlet_tube_ID,
                              mdot = mdot_guess, 
                              State1 = inletState.copy(),
                              fixed = 1,
                              TubeFcn = sim.TubeCode))
            sim.add_tube(Tube(key1 = 'outlet.1',
                              key2 = 'outlet.2',
                              L = sim.outlet_tube_length,
                              ID = sim.outlet_tube_ID,
                              mdot = mdot_guess, 
                              State2 = outletState.copy(),
                              fixed = 2,
                              TubeFcn = sim.TubeCode))
                                     
            # Add all the leakage flows
            sim.auto_add_leakage(flankFunc = sim.FlankLeakage, 
                                 radialFunc = sim.RadialLeakage)
                                 
            # Add the inlet-to-shell flow with a fixed area
            FP = FlowPath(key1='inlet.2',
                  key2='sa', 
                  MdotFcn=IsentropicNozzleWrapper(),
                  )
            FP.A = pi*sim.inlet_tube_ID**2/4*{Xd_inlet:s}
            sim.add_flow(FP)
            
            # Add the suction-area to suction chambers flows
            sim.add_flow(FlowPath(key1='sa', 
                                  key2='s1',
                                  MdotFcn=sim.SA_S1,
                                  MdotFcn_kwargs = dict(X_d = {Xd_sa_s1:s})
                                  )
                        )
            sim.add_flow(FlowPath(key1 = 'sa',
                                  key2 = 's2',
                                  MdotFcn = sim.SA_S2,
                                  MdotFcn_kwargs = dict(X_d = {Xd_sa_s2:s})
                                  )
                        )
            
            FP = FlowPath(key1='outlet.1',
                          key2='dd',
                          MdotFcn=IsentropicNozzleWrapper(),
                          )
            FP.A = pi*sim.d_discharge**2/4
            sim.add_flow(FP)
            
            FP = FlowPath(key1='outlet.1', 
                          key2='ddd', 
                          MdotFcn=IsentropicNozzleWrapper(),
                          )
            FP.A = pi*sim.d_discharge**2/4
            sim.add_flow(FP)
            
            sim.add_flow(FlowPath(key1='d1',
                                  key2='dd',
                                  MdotFcn=sim.D_to_DD))
            sim.add_flow(FlowPath(key1='d2',
                                  key2='dd',
                                  MdotFcn=sim.D_to_DD))
            """.format(Xd_sa_s1 = Xd_dict['sa-s1'],
                       Xd_sa_s2 = Xd_dict['sa-s2'],
                       Xd_inlet = Xd_dict['inlet.2-sa'],
                       )
            )
        
       
#    def collect_output_terms(self):
#        _T = []
#        
#        for i,Tube in zip(:
#            _T.extend([dict(attr = "Tubes["+str(i)+"].State1.T",
#                            text = "Tube T ["+ str(Tube.key1) +"] [K]",
#                            parent = self
#                            ),
#                       dict(attr = "Tubes["+str(i)+"].State2.T",
#                            text = "Tube T ["+ str(Tube.key2) +"] [K]",
#                            parent = self
#                            ),
#                       dict(attr = "Tubes["+str(i)+"].State1.p",
#                            text = "Tube p ["+ str(Tube.key1) +"] [kPa]",
#                            parent = self
#                            ),
#                       dict(attr = "Tubes["+str(i)+"].State2.p",
#                            text = "Tube p ["+ str(Tube.key2) +"] [kPa]",
#                            parent = self
#                            )
#                       ])
#        return _T
    
class MechanicalLossesPanel(pdsim_panels.PDPanel):
    
    desc_map = dict(h_shell = ('Shell-ambient mean HTC','W/m^2/K'),
                    A_shell = ('Shell outer area [m\xb2]','m^2'),
                    Tamb = ('Ambient temperature [K]','K'),
                    
                    mu_oil = ('Viscosity of the oil [Pa-s]','Pa-s'),
                    D_upper_bearing = ('Upper bearing journal diameter [m]','m'),
                    L_upper_bearing = ('Upper bearing length [m]','m'),
                    c_upper_bearing = ('Upper bearing clearance [m]','m'),
                    D_crank_bearing = ('Crank bearing journal diameter [m]','m'),
                    L_crank_bearing = ('Crank bearing length [m]','m'),
                    c_crank_bearing = ('Upper bearing clearance [m]','m'),
                    D_lower_bearing = ('Lower bearing journal diameter [m]','m'),
                    L_lower_bearing = ('Lower bearing length [m]','m'),
                    c_lower_bearing = ('Upper bearing clearance [m]','m'),
                    journal_tune_factor = ('Tuning factor on journal bearing losses [-]','-'),
                    thrust_friction_coefficient = ('Thrust bearing friction coefficient [-]','-'),
                    thrust_ID = ('Thrust bearing inner diameter [m]','m'),
                    thrust_OD = ('Thrust bearing outer diameter [m]','m'), 
                    orbiting_scroll_mass = ('Orbiting scroll mass [kg]','kg'), 
                    L_ratio_bearings = ('Ratio of lengths to the bearings [-]','-'),
                    
                    HTC = ('Heat transfer coefficient in the scrolls [W/m\xb2/K]','W/m^2/K'),
                    )
    
    def __init__(self, parent, config, **kwargs):
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        # Now we are going to put everything into a scrolled window
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # The scrolled panel
        scrolled_panel = ScrolledPanel(self, size = (-1,-1), style = wx.TAB_TRAVERSAL, name="panel1")
        scrolled_panel.SetScrollbars(1, 1, 1, 1)
        
        # The sizer for all the objects
        sizer_for_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        """
        There are 2 possibilities for the types of motor models supported.
        
        The motor can be map based in which case efficiency and slip speed are
        given as a function of mechanical torque output.  Or the efficiency and
        rotational speed are given
         
        Either the motor rejects its heat to the ambient (as in open-drive), or
        it rejects its heat to the suction volume
        """
        
        self.motor_choices = MotorChoices(scrolled_panel)
        
        if ('eta_motor' in config
            and 'eta_motor_coeffs' not in config
            and 'tau_motor_coeffs' not in config
            and 'omega_motor_coeffs' not in config):
            #Only eta_motor is provided, use it in the motor panel
            self.motor_choices.SetSelection(0)
            #Set the value in the panel
            self.motor_choices.eta_motor.SetValue(str(config['eta_motor']))
            
        elif ('eta_motor' not in config
            and 'eta_motor_coeffs' in config
            and 'tau_motor_coeffs' in config
            and 'omega_motor_coeffs' in config):
            #Coefficients are provided, use them in the motor panel
            self.motor_choices.SetSelection(1)
            values = [config['tau_motor_coeffs'],
                      config['eta_motor_coeffs'],
                      config['omega_motor_coeffs']
                      ]
            self.motor_choices.MCT.update_from_configfile(values)
        else:
            raise ValueError('Your combination of motor terms is not valid')
        
        keys = ['h_shell','A_shell','Tamb','mu_oil',
                'D_upper_bearing','L_upper_bearing','c_upper_bearing',
                'D_crank_bearing','L_crank_bearing','c_crank_bearing',
                'D_lower_bearing','L_lower_bearing','c_lower_bearing',
                'journal_tune_factor',
                'thrust_friction_coefficient', 'thrust_ID', 'thrust_OD', 
                'orbiting_scroll_mass', 'L_ratio_bearings', 'HTC'
                ]
        
        # The list for all the annotated objects
        self.annotated_values = []
        
        # Loop over the first group of inputs
        for key in keys:
            # Get the annotation and the units for the term 
            annotation, units = self.desc_map[key]
            # Add the annotated object to the list of objects
            self.annotated_values.append(AnnotatedValue(key, config[key], annotation, units))
            
        # Build the items and return the list of annotated GUI objects
        annotated_GUI_objects = self.construct_items(self.annotated_values, 
                                                     sizer_for_inputs, 
                                                     parent = scrolled_panel)
        
        # Register terms in the GUI database
        self.main.register_GUI_objects(annotated_GUI_objects)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(wx.StaticText(scrolled_panel, -1, "Motor Model"))
        sizer.Add(wx.StaticLine(scrolled_panel, -1, (25, 50), (300,1)))
        sizer.Add(self.motor_choices)
        sizer.AddSpacer(20)
        sizer.Add(sizer_for_inputs)
        
        scrolled_panel.SetSizer(sizer)
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND)
        self.SetSizer(main_sizer)
        
        sizer.Layout()
        
    def get_script_chunks(self):
        """
        Returns a formatted string for the script that will be execfile-d
        """
        if self.motor_choices.GetSelection() == 0:
            #Use the value for the motor efficiency
            return textwrap.dedent(
               """
               sim.motor = Motor()
               sim.motor.set_eta({eta_motor:s})
               sim.motor.suction_fraction = 1.0 #hard-coded
               """.format(eta_motor = self.motor_choices.eta_motor.GetValue()))
        elif self.motor_choices.GetSelection() == 1:
            # Get the tuple of list of coeffs from the MCT, then unpack the tuple
            # back into the call to set the coefficients
            c = self.motor_choices.MCT.get_coeffs()
            #Will set the type flag itself
            return textwrap.dedent(
               """
               sim.motor = Motor()
               sim.motor.set_coeffs(tau_coeffs = {tau_coeffs:s},
                                    eta_coeffs = {eta_coeffs:s},
                                    omega_coeffs = {omega_coeffs:s})
               sim.motor.suction_fraction = 1.0 #hard-coded
               """.format(tau_coeffs = str(c[0]),
                          eta_coeffs = str(c[1]),
                          omega_coeffs = str(c[2])
                          )        )
        else:
            raise NotImplementedError
        
    def post_get_from_configfile(self,key,value):
        """
        Parse the coefficients or other things that are not handled in the normal way
        """
        #value starts as a string like 'coeffs,1;2;3' --> list of string
        values = value.split(',')[1].split(';')
        #Return a list of coefficients - whitespace is ignored in the float() call
        return [float(v) for v in values] 
        
    def collect_output_terms(self):
        return [dict(attr = "losses.crank_bearing",
                       text = "Crank bearing loss [kW]",
                       parent = self
                       ),
                dict(attr = "losses.upper_bearing",
                       text = "Upper bearing loss [kW]",
                       parent = self
                       ),
                dict(attr = "losses.lower_bearing",
                       text = "Lower bearing loss [kW]",
                       parent = self
                       ),
                dict(attr = "losses.thrust_bearing",
                       text = "Thrust bearing loss [kW]",
                       parent = self
                       ),
                dict(attr = 'forces.mean_Fm',
                     text = 'Mean pin force magnitude [kN]',
                     parent = self),
                dict(attr = 'forces.mean_Fz',
                     text = 'Mean axial force [kN]',
                     parent = self),
                dict(attr = 'forces.mean_Fr',
                     text = 'Mean radial force [kN]',
                     parent = self),
                dict(attr = 'forces.mean_Ft',
                     text = 'Mean tangential force [kN]',
                     parent = self),
                dict(attr = 'forces.inertial',
                     text = 'Inertial forces on orb. scroll [kN]',
                     parent = self),
                  ]
    
    def post_prep_for_configfile(self):
        """
        Custom code for output to config file to handle motor
        """
        if self.motor_choices.GetSelection() == 0:
            eta = self.motor_choices.eta_motor.GetValue()
            return 'eta_motor = float, '+eta+', Motor Efficiency [-]'
        else:
            return self.motor_choices.MCT.string_for_configfile()
#            
#        
#if __name__=='__main__':
#    app = wx.App(False)
#    
#    frame = wx.Frame()
#    frame.geo = GeometryPanel(frame,) 
#    frame.Show(True) 
#    
#    app.MainLoop()
        