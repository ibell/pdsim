# -*- coding: latin-1 -*-
import pdsim_panels
import wx
from math import pi
import textwrap

from datatypes import HeaderStaticText
from wx.lib.scrolledpanel import ScrolledPanel

class GeometryPanel(pdsim_panels.PDPanel):
    
    desc_map = dict(piston_diameter = ('Piston Diameter [m]','m'),
                    piston_length = ('Piston Length [m]','m'),
                    crank_length = ('Crank length [m]','m'),
                    connecting_rod_length = ('Connecting rod length [m]','m'),
                    x_TDC = ('Distance to piston at TDC [m]','m'),
                    shell_volume = ('Shell volume [m\xb3]','m^3'),
                    inlet_tube_length = ('Inlet tube length [m]','m',0.02),
                    inlet_tube_ID = ('Inlet tube inner diameter [m]','m',0.02),
                    outlet_tube_length = ('Outlet tube length [m]','m',0.02),
                    outlet_tube_ID = ('Outlet tube inner diameter [m]','m',0.02),
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
        keys = ['piston_diameter','piston_length','crank_length',
                'connecting_rod_length','x_TDC','shell_volume',
                'inlet_tube_length','inlet_tube_ID','outlet_tube_length',
                'outlet_tube_ID']
        annotated_values = self.get_annotated_values(keys)
            
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
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND)#|wx.ALIGN_CENTER_HORIZONTAL)
        self.SetSizer(main_sizer)
        
        sizer.Layout()
        
    def get_script_chunks(self, plugin_chunks = None):
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
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND)#|wx.ALIGN_CENTER_HORIZONTAL)
        self.SetSizer(main_sizer)
        
        sizer.Layout()
        
    def get_script_chunks(self, plugin_chunks = None):
        chunk = ''
        
        def get(key):
            # Compact code to get a parameter from the main database
            return self.main.get_GUI_object(key).GetValue()
        
        keys = ['valve_d','valve_C_D','valve_h','valve_a','valve_l','valve_rho',
                'valve_x_stopper','valve_E', 'inlet_tube_length','inlet_tube_ID',
                'outlet_tube_length','outlet_tube_ID']
        
        for term in keys:
            val = self.main.get_GUI_object_value(term)
            chunk += 'sim.{name:s} = {value:s}\n'.format(name = term, value = str(val))        
        
        for term in ['d_discharge','d_suction']:
            val = self.main.get_GUI_object_value(term)
            chunk += 'sim.{name:s} = {value:s}\n'.format(name = term, value = str(val))

        chunk += textwrap.dedent(
        """
        
        # First add the control volumes.
        sim.add_CV( ControlVolume(key = 'A',
                                  initialState = outletState.copy(),
                                  VdVFcn = sim.V_dV,
                                  becomes = 'A') )
        sim.add_CV( ControlVolume(key = 'shell',
                                  initialState = inletState.copy(),
                                  VdVFcn = sim.V_shell,
                                  becomes = 'shell') )
                                  
        sim.add_flow(FlowPath(key1='shell', key2='inlet.2', MdotFcn = sim.Inlet))
        sim.add_flow(FlowPath(key1='shell', key2='A', MdotFcn = sim.Suction))
        sim.add_flow(FlowPath(key1='outlet.1', key2='A', MdotFcn = sim.Discharge))
        sim.add_flow(FlowPath(key1='shell', key2='A', MdotFcn = sim.PistonLeakage))
        
        # Calculate Vdisp
        sim.pre_solve()
        sim.Vdisp = sim.Vdisp()
    
        # Get the guess for the mass flow rate
        mdot_guess = inletState.rho*sim.Vdisp*sim.omega/(2*pi)

        # Add both the inlet and outlet tubes
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
                              
        E = sim.valve_E
        I = (sim.valve_d*sim.valve_h**3)/12                   #Moment of Intertia for valve,[m^4]
        k_valve = (6*E*I)/(sim.valve_a**2*(3*sim.valve_l-sim.valve_a))   #Valve stiffness
        m_eff = (1/3)*sim.valve_rho*sim.valve_l*sim.valve_d*sim.valve_h      #Effective mass of valve reeds
        x_tr_suction = 0.25*(sim.d_suction**2/sim.valve_d)
        x_tr_discharge = 0.25*(sim.d_discharge**2/sim.valve_d)
        
        #  The suction valve parameters
        sim.suction_valve = ValveModel(
              d_valve = sim.valve_d,
              d_port = sim.d_suction,
              C_D = sim.valve_C_D,
              rho_valve = sim.valve_rho,
              x_stopper = sim.valve_x_stopper,
              m_eff = m_eff,
              k_valve = k_valve,
              x_tr = x_tr_suction,
              key_up = 'inlet.2',
              key_down = 'A'
              )
        sim.add_valve(sim.suction_valve)
        
        #  The discharge valve parameters
        sim.discharge_valve=ValveModel(
              d_valve = sim.valve_d,
              d_port = sim.d_discharge,
              C_D = sim.valve_C_D,
              rho_valve = sim.valve_rho,
              x_stopper = sim.valve_x_stopper,
              m_eff = m_eff,
              k_valve = k_valve,
              x_tr = x_tr_discharge,
              key_up='A',
              key_down='outlet.1'
              )
        sim.add_valve(sim.discharge_valve)
        """).format()
        
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
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND)#|wx.ALIGN_CENTER_HORIZONTAL)
        self.SetSizer(main_sizer)
        
        sizer.Layout()
        
    def get_script_chunks(self, plugin_chunks = None):
        chunk = ''
        for term in ['eta_motor','h_shell','A_shell','Tamb','mu_oil','delta_gap']:
            val = self.main.get_GUI_object_value(term)
            chunk += 'sim.{name:s} = {value:s}\n'.format(name = term, value = str(val))
        return chunk