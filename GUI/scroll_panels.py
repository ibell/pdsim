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
from CoolProp import State as CPState

LabeledItem = pdsim_panels.LabeledItem

class PlotPanel(wx.Panel):
    def __init__(self, parent, **kwargs):
        size = kwargs.get('size',(200,200))
        wx.Panel.__init__(self, parent, size = size, **kwargs)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.figure = mpl.figure.Figure(dpi=100, figsize=(size[0]/100, size[1]/100))
        self.canvas = WXCanvas(self, -1, self.figure)
#        self.canvas.resize(200,200)
        #self.toolbar = WXToolbar(self.canvas)
        #self.toolbar.Realize()
        sizer.Add(self.canvas)
        #sizer.Add(self.toolbar)
        self.SetSizer(sizer)
        sizer.Layout()
        

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
        SAF = ScrollAnimForm(self.Scroll.geo, size=(400,400))
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
        self.Scroll.set_disc_geo('2Arc', r2=0)
        
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
        phi_i0 = float(self.phi_fi0.GetValue())
        phi_is = float(self.phi_fis.GetValue())
        phi_os = float(self.phi_fos.GetValue())
        
        #Set the scroll wrap geometry
        scroll.set_scroll_geo(Vdisp,Vratio,t,ro,
                                   phi_i0 = phi_i0, 
                                   phi_os = phi_os,
                                   phi_is = phi_is)
        scroll.set_disc_geo('2Arc', r2 = 0)
        
        # self.delta_flank and self.delta_radial are set in the conventional way
        # using the parametric table
        scroll.geo.delta_flank = scroll.delta_flank
        scroll.geo.delta_radial = scroll.delta_radial
        
        if self.UseOffset.IsChecked():
            scroll.geo.phi_ie_offset = pi
        else:
            scroll.geo.phi_ie_offset = 0
            
        print self.post_set_params_str()
        
    def post_set_params_str(self):
    
        #Parameters to be set in the string:
        str_params = dict(Vdisp = self.Vdisp.GetValue(),
                          Vratio = self.Vratio.GetValue(),
                          t = self.t.GetValue(),
                          ro = self.ro.GetValue(),
                          phi_i0 = self.phi_fi0.GetValue(),
                          phi_is = self.phi_fis.GetValue(),
                          phi_os = self.phi_fos.GetValue())

        return textwrap.dedent(
            """
            #Parameters from the GUI
            Vdisp = {Vdisp:s} #[m]
            Vratio = {Vratio:s} #[-] 
            t = {t:s} #[m]
            ro = {ro:s} #[m]
            phi_i0 = {phi_i0:s} #[rad]
            phi_is = {phi_is:s} #[rad]
            phi_os = {phi_os:s} #[rad]
            
            #Set the scroll wrap geometry
            scroll.set_scroll_geo(Vdisp,Vratio,t,ro,
                                       phi_i0 = phi_i0, 
                                       phi_os = phi_os,
                                       phi_is = phi_is)
            scroll.set_disc_geo('2Arc', r2 = 0)
            
            # self.delta_flank and self.delta_radial are set in the conventional way
            # using the parametric table
            scroll.geo.delta_flank = scroll.delta_flank
            scroll.geo.delta_radial = scroll.delta_radial
            
            if self.UseOffset.IsChecked():
                scroll.geo.phi_ie_offset = pi
            else:
                scroll.geo.phi_ie_offset = 0        
            """.format(**str_params))
    
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
#        Xd_sa_s1 = self.get_from_configfile('MassFlowPanel','Xd_sa_s1')
#        self.suctionflow1.set_attr('X_d', float(Xd_sa_s1))
        
        self.suctionflow2 = SuctionFlowChoice(parent = self,
                                              key1 = 'sa',
                                              key2 = 's2',
                                              label = 'Flow to suction chamber #2')    
#        Xd_sa_s2 = self.get_from_configfile('MassFlowPanel','Xd_sa_s2')
#        self.suctionflow2.set_attr('X_d', float(Xd_sa_s2))
        
        self.inletflow = InletFlowChoice(parent = self,
                                         key1 = 'inlet.2',
                                         key2 = 'sa',
                                         label = 'Flow into shell',
                                         )
#        Xd_shell_sa = self.get_from_configfile('MassFlowPanel','Xd_shell_sa')
#        self.inletflow.set_attr('X_d', float(Xd_shell_sa))
        
        self.flows = [self.suctionflow1, self.suctionflow2, self.inletflow]
        
        box_sizer.AddMany(self.flows)
        
        self.resize_flows(self.flows)
        
        self.SetSizer(box_sizer)
        box_sizer.Layout()
        
        self.items=self.items1
        
    def resize_flows(self, flows):
        """
        Resize the labels for the flows to all be the same size
        """
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
                A = pi * D**2/4 * param_dict.pop('X_d')
                param_dict['A']=A
                
            simulation.add_flow(FlowPath(key1 = flow.key1,
                                         key2 = flow.key2,
                                         MdotFcn = func,
                                         MdotFcn_kwargs = param_dict
                                         )
                                )
            
        A_discharge_port = pi*simulation.d_discharge**2/4
        
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
                                  mdot=simulation.inletState.copy().rho*Vdot, 
                                  State2=outletState.copy(),
                                  fixed=2,
                                  TubeFcn=simulation.TubeCode) )
        
       
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
     
class MotorCoeffsTable(wx.ListCtrl, TextEditMixin):
    
    def __init__(self, parent, values = None):
        """
        Parameters
        ----------
        parent : wx.window
            The parent of this checklist
        values : A 3-element list of lists for all the coeff
        """
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT)
        TextEditMixin.__init__(self)
        
        #Build the headers
        self.InsertColumn(0,'Motor Torque [N-m]')
        self.InsertColumn(1,'Efficiency [-]')
        self.InsertColumn(2,'Slip speed [rad/s]')
        
        #: The values stored as a list of lists in floating form
        self.values = values
        
        #Reset the values
        self.refresh_table()
        
        #Set the column widths    
        for i in range(3):
            self.SetColumnWidth(i,wx.LIST_AUTOSIZE_USEHEADER)
        
        # Turn on callback to write values back into internal data structure when
        # a cell is edited
        self.Bind(wx.EVT_LIST_END_LABEL_EDIT, self.OnCellEdited)
        
        #Required width for the table
        min_width = sum([self.GetColumnWidth(i) for i in range(self.GetColumnCount())])
        
        #No required height (+30 for buffer to account for vertical scroll bar)
        self.SetMinSize((min_width + 30,-1))

    def OnCellEdited(self, event):
        """
        Once the cell is edited, write its value back into the data matrix
        """
        row_index = event.m_itemIndex
        col_index = event.Column
        val = float(event.Text)
        self.data[row_index][col_index-1] = val
    
    def GetStringCell(self,Irow,Icol):
        """ Returns a string representation of the cell """
        return self.data[Irow][Icol]
    
    def GetFloatCell(self,Irow,Icol):
        """ Returns a float representation of the cell """
        return float(self.data[Irow][Icol])
    
    def AddRow(self):
        
        row = [0]*self.GetColumnCount()
        
        i = len(self.data)-1
        self.InsertStringItem(i,'')
        for j,val in enumerate(row):
            self.SetStringItem(i,j+1,str(val))
        self.CheckItem(i)
        
        self.data.append(row)
        
    def RemoveRow(self, i = 0):
        self.data.pop(i)
        self.DeleteItem(i)
        
    def RemoveLastRow(self):
        i = len(self.data)-1
        self.data.pop(i)
        self.DeleteItem(i)
    
    def update_from_configfile(self, values):
        """
        
        Parameters
        ----------
        values : list of lists, with entries as floating point values
            The first entry is a list (or other iterable) of torque values
            
            The second entry is a list (or other iterable) of efficiency values
            
            The third entry is a list (or other iterable) of slip speed values
        """
        self.values = values
        self.refresh_table()
        
    def string_for_configfile(self):
        """
        Build and return a string for writing to the config file
        """
            
        tau_list = self.values[0]
        tau_string = 'tau_motor_coeffs = coeffs, '+'; '.join([str(tau) for tau in tau_list])
        
        eta_list = self.values[1]
        eta_string = 'eta_motor_coeffs = coeffs, '+'; '.join([str(eta) for eta in eta_list])
        
        omega_list = self.values[2]
        omega_string = 'omega_motor_coeffs = coeffs, '+'; '.join([str(omega) for omega in omega_list])
            
        return tau_string + '\n' + eta_string + '\n' + omega_string + '\n'
        
        
    def refresh_table(self):
        """
        Take the values from self.values and write them to the table
        """
        #Remove all the values in the table
        for i in reversed(range(self.GetItemCount())):
            self.DeleteItem(i)
            
        if self.values is None:
            #Add a few rows
            for i in range(10):
                self.InsertStringItem(i,str(i))
        else:
            #They all need to be the same length
            assert len(self.values[0]) == len(self.values[1]) == len(self.values[2])
            for i in range(len(self.values[0])):
                self.InsertStringItem(i,str(self.values[0][i]))
                self.SetStringItem(i,1,str(self.values[1][i]))
                self.SetStringItem(i,2,str(self.values[2][i]))
                
    def get_coeffs(self):
        """
        Get the list of lists of values that are used in the table
        """
        return self.values
    

        
class MotorChoices(wx.Choicebook):
    def __init__(self, parent):
        wx.Choicebook.__init__(self, parent, -1)
        
        self.pageconsteta=wx.Panel(self)
        self.AddPage(self.pageconsteta,'Constant efficiency')
        self.eta_motor_label, self.eta_motor = LabeledItem(self.pageconsteta, 
                                                           label="Motor Efficiency [-]",
                                                           value='0.9')
        sizer=wx.FlexGridSizer(cols = 2, hgap = 3, vgap = 3)
        sizer.AddMany([self.eta_motor_label, self.eta_motor])
        self.pageconsteta.SetSizer(sizer)
        
        self.pagemotormap=wx.Panel(self)
        self.AddPage(self.pagemotormap,'Motor map')
        self.MCT = MotorCoeffsTable(self.pagemotormap,values = [[1,2,3],[0.9,0.9,0.9],[307,307,307]])
        sizer=wx.FlexGridSizer(cols = 2, hgap = 3, vgap = 3)
        sizer.Add(self.MCT, 1, wx.EXPAND)
        self.pagemotormap.SetSizer(sizer)
        sizer.Layout()
    
class MechanicalLossesPanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, configfile,**kwargs):
    
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        #Loads all the parameters from the config file (case-sensitive)
        self.configdict, self.descdict = self.get_from_configfile('MechanicalLossesPanel')
        
        #: The box sizer that contains all the sizers
        box_sizer = wx.BoxSizer(wx.VERTICAL)
        
        """
        There are 2 possibilities for the types of motor models supported.
        
        The motor can be map based in which case efficiency and slip speed are
        given as a function of mechanical torque output.  Or the efficiency and
        rotational speed are given
         
        Either the motor rejects its heat to the ambient (as in open-drive), or
        it rejects its heat to the suction volume
        """
        
        box_sizer.Add(wx.StaticText(self, -1, "Motor Model"))
        box_sizer.Add(wx.StaticLine(self, -1, (25, 50), (300,1)))
        
        self.motor_choices = MotorChoices(self)
        box_sizer.Add(self.motor_choices)
        
        if ('eta_motor' in self.configdict 
            and 'eta_motor_coeffs' not in self.configdict
            and 'tau_motor_coeffs' not in self.configdict
            and 'omega_motor_coeffs' not in self.configdict):
            #Only eta_motor is provided, use it in the motor panel
            self.motor_choices.SetSelection(0)
            #Set the value in the panel
            self.motor_choices.eta_motor.SetValue(str(self.configdict['eta_motor']))
            
        elif ('eta_motor' not in self.configdict 
            and 'eta_motor_coeffs' in self.configdict
            and 'tau_motor_coeffs' in self.configdict
            and 'omega_motor_coeffs' in self.configdict):
            #Coefficients are provided, use them in the motor panel
            self.motor_choices.SetSelection(1)
            values = [self.configdict['tau_motor_coeffs'],
                      self.configdict['eta_motor_coeffs'],
                      self.configdict['omega_motor_coeffs']
                      ]
            self.motor_choices.MCT.update_from_configfile(values)
        elif False:
            #If neither of the above are matched, fall back to the default configuration
            pass
        else:
            raise ValueError('Your combination of motor terms is not valid')
        
        self.items = [
      
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
                    dict(attr='journal_tune_factor'),
                    
                    dict(attr = 'thrust_friction_coefficient'),
                    dict(attr = 'thrust_ID'),
                    dict(attr = 'thrust_OD'),
                    dict(attr = 'orbiting_scroll_mass'),
                    dict(attr = 'L_ratio_bearings'),
                    dict(attr = 'HTC')
        ]
        
        
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
#        self.ConstructItems([self.items[0]],sizer,self.configdict,self.descdict)
#        
#        box_sizer.AddSpacer(20)
        self.ConstructItems(self.items,sizer,self.configdict,self.descdict)
        
        box_sizer.Add(sizer)
        self.SetSizer(box_sizer)
        sizer.Layout()
        
    def post_set_params(self, sim):
        """
        Set other parameters after loading all the conventional terms
        
        Parameters
        ----------
        sim : PDSimCore instance (or derived class)
            The simulation instance
        """
        #Create the motor if if doesn't have one
        if not hasattr(sim,'motor'):
            sim.motor = Motor()
        
        if self.motor_choices.GetSelection() == 0:
            #Use the value for the motor efficiency
            sim.eta_motor = float(self.motor_choices.eta_motor.GetValue())
            sim.motor.set_eta(sim.eta_motor)
        elif self.motor_choices.GetSelection() == 1:
            # Get the tuple of list of coeffs from the MCT, then unpack the tuple
            # back into the call to set the coefficients
            c = self.motor_choices.MCT.get_coeffs()
            #Set the coefficients for the motor model
            sim.motor.set_coeffs(tau_coeffs = c[0], eta_coeffs = c[1], omega_coeffs = c[2])
            #Will set the type flag itself
        else:
            raise NotImplementedError
        
        sim.motor.suction_fraction = 1.0
        
            
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
            
        
        