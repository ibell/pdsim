# -*- coding: latin-1 -*-
import pdsim_panels
import wx
from math import pi

from datatypes import HeaderStaticText
from wx.lib.scrolledpanel import ScrolledPanel

class GeometryPanel(pdsim_panels.PDPanel):
    
    desc_map = dict(piston_diameter = ('Piston Diameter [m]','m'),
                    piston_length = ('Piston Length [m]','m'),
                    crank_length = ('Crank length [m]','m'),
                    connecting_rod_length = ('Connecting rod length [m]','m'),
                    x_TDC = ('Distance to piston at TDC [m]','m'),
                    shell_volume = ('Shell volume [m\xb3]','m^3')
                    )
    
    def __init__(self, parent, config, **kwargs):
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        # Now we are going to put everything into a scrolled window
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # The scrolled panel
        scrolled_panel = ScrolledPanel(self, size = (-1,-1), style = wx.TAB_TRAVERSAL, name="panel1")
        scrolled_panel.SetScrollbars(1, 1, 1, 1)
        
        annotated_GUI_objects = []
        self.config = config
        self.keys_for_config = []
        
        #----------------------------------------------------------------------
        # The sizer for all the heat transfer terms
        sizer_for_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the HT inputs
        annotated_values = self.get_annotated_values(['piston_diameter','piston_length','crank_length','connecting_rod_length','x_TDC','shell_volume'])
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                       sizer = sizer_for_inputs,
                                                       parent = scrolled_panel)
        
        # Register terms in the GUI database
        self.main.register_GUI_objects(annotated_GUI_objects)
        
        #self.main.get_GUI_object('L_ratio_bearings').GUI_location.SetToolTipString('Ratio of z1/z2, where\n\nz1 : the length from the centerline of the upper bearing to the lower bearing\nz2 : the length from the centerline of the upper bearing to the orbiting scroll bearing')
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(HeaderStaticText(scrolled_panel, "Geometric Inputs"), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(sizer_for_inputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        scrolled_panel.SetSizer(sizer)
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)
        self.SetSizer(main_sizer)
        
        sizer.Layout()
        
    def get_script_chunks(self):
        chunk = ''
        for term in ['piston_diameter','piston_length','crank_length','connecting_rod_length','x_TDC','shell_volume']:
            val = self.main.get_GUI_object_value(term)
            chunk += 'sim.{name:s} = {value:s}\n'.format(name = term, value = str(val))
        return chunk
        
class MassFlowPanel(pdsim_panels.PDPanel):
    
    desc_map = dict(d_discharge = ('Discharge port diameter [m]','m',0.0059),
                    d_suction = ('Suction port diameter [m]','m',0.0059),
                    valve_E = ('Youngs Modulus [Pa]','Pa',1.93e+11),
                    valve_d = ('Valve diameter [m]','m',0.007),
                    valve_h = ('Valve thickness [m]','m',0.0001532),
                    valve_l = ('Valve length [m]','m',0.018),
                    valve_a = ('Valve distance from anchor [m]','m',0.014),
                    valve_x_stopper = ('Valve distance to stopper [m]','m',0.0018),
                    valve_rho = ('Valve metal density [kg/m\xb3]','kg/m^3',8000.0),
                    valve_C_D = ('Valve drag coefficient [-]','-',1.17),
                    )
    
    def __init__(self, parent, config, **kwargs):
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        # Now we are going to put everything into a scrolled window
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # The scrolled panel
        scrolled_panel = ScrolledPanel(self, size = (-1,-1), style = wx.TAB_TRAVERSAL, name="panel1")
        scrolled_panel.SetScrollbars(1, 1, 1, 1)
        
        annotated_GUI_objects = []
        self.config = config
        self.keys_for_config = []
        
        #----------------------------------------------------------------------
        # The sizer for all the heat transfer terms
        sizer_for_port_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the HT inputs
        annotated_values = self.get_annotated_values(['d_discharge','d_suction'])
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                       sizer = sizer_for_port_inputs,
                                                       parent = scrolled_panel)
                                                       
       #----------------------------------------------------------------------
        # The sizer for all the valve terms
        sizer_for_valve_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
    
        # Loop over the HT inputs
        annotated_values = self.get_annotated_values(['valve_E','valve_d','valve_h','valve_l','valve_a','valve_x_stopper','valve_rho','valve_C_D'])
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                       sizer = sizer_for_valve_inputs,
                                                       parent = scrolled_panel)
                                                       
                                                       
        
        # Register terms in the GUI database
        self.main.register_GUI_objects(annotated_GUI_objects)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(HeaderStaticText(scrolled_panel, "Geometric Inputs"), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(sizer_for_port_inputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        sizer.Add(HeaderStaticText(scrolled_panel, "Valve Inputs"), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(sizer_for_valve_inputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        scrolled_panel.SetSizer(sizer)
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)
        self.SetSizer(main_sizer)
        
        sizer.Layout()
        
    def get_script_chunks(self):
        chunk = ''
        for term in ['d_discharge','d_suction','valve_E','valve_d',
                     'valve_h','valve_l','valve_a','valve_x_stopper','valve_rho',
                     'valve_C_D']:
            val = self.main.get_GUI_object_value(term)
            chunk += 'sim.{name:s} = {value:s}\n'.format(name = term, value = str(val))
        return chunk
        
class MechanicalLossesPanel(pdsim_panels.PDPanel):
    
    desc_map = dict(eta_motor = ('Motor efficiency [-]','-',0.95),
                    h_shell = ('Shell air-side heat transfer coefficient [kW/m\xb2/K]','kW/m^2/K',0.01),
                    A_shell = ('Shell Area [m\xb2]','m^2','0.040536'),
                    Tamb = ('Ambient temperature [K]','K',298.0),
                    mu_oil = ('Oil viscosity [Pa-s]','Pa-s',0.0086),
                    delta_gap = ('Gap width [m]','m',2e-5), 
                    )
    
    def __init__(self, parent, config, **kwargs):
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        # Now we are going to put everything into a scrolled window
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # The scrolled panel
        scrolled_panel = ScrolledPanel(self, size = (-1,-1), style = wx.TAB_TRAVERSAL, name="panel1")
        scrolled_panel.SetScrollbars(1, 1, 1, 1)
        
        annotated_GUI_objects = []
        self.config = config
        self.keys_for_config = []
        
        #----------------------------------------------------------------------
        # The sizer for all the heat transfer terms
        sizer_for_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the inputs
        annotated_values = self.get_annotated_values(['eta_motor','h_shell','A_shell','Tamb','mu_oil','delta_gap'])
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                       sizer = sizer_for_inputs,
                                                       parent = scrolled_panel)
                                                       
        
        # Register terms in the GUI database
        self.main.register_GUI_objects(annotated_GUI_objects)
        
        #self.main.get_GUI_object('L_ratio_bearings').GUI_location.SetToolTipString('Ratio of z1/z2, where\n\nz1 : the length from the centerline of the upper bearing to the lower bearing\nz2 : the length from the centerline of the upper bearing to the orbiting scroll bearing')
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(HeaderStaticText(scrolled_panel, "Geometric Inputs"), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(sizer_for_inputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        scrolled_panel.SetSizer(sizer)
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)
        self.SetSizer(main_sizer)
        
        sizer.Layout()
        
    def get_script_chunks(self):
        chunk = ''
        for term in ['eta_motor','h_shell','A_shell','Tamb','mu_oil','delta_gap']:
            val = self.main.get_GUI_object_value(term)
            chunk += 'sim.{name:s} = {value:s}\n'.format(name = term, value = str(val))
        return chunk