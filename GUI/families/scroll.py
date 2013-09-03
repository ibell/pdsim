from __future__ import absolute_import
import wx, yaml, os
from panels import scroll_panels, pdsim_panels

family_menu_name = 'Scroll Compressor'

class InputsToolBook(pdsim_panels.InputsToolBook):
    """
    The toolbook that contains the pages with input values
    """
    def __init__(self, parent, config):
        """
        Parameters
        ----------
        config : yaml config file
            The top-level configuration file
        """
        wx.Toolbook.__init__(self, parent, -1, style=wx.BK_LEFT)
        il = wx.ImageList(32, 32)
        indices=[]
        for imgfile in ['Geometry.png',
                        'StatePoint.png',
                        'MassFlow.png',
                        'MechanicalLosses.png']:
            ico_path = os.path.join('ico',imgfile)
            indices.append(il.Add(wx.Image(ico_path,wx.BITMAP_TYPE_PNG).ConvertToBitmap()))
        self.AssignImageList(il)
        
        Main = wx.GetTopLevelParent(self)
        # Make the scroll panels.  
        self.panels=(scroll_panels.GeometryPanel(self, config['GeometryPanel'],name='GeometryPanel'),
                     pdsim_panels.StateInputsPanel(self, config['StatePanel'], name='StatePanel'),
                     scroll_panels.MassFlowPanel(self, config['MassFlowPanel'], name='MassFlowPanel'),
                     scroll_panels.MechanicalLossesPanel(self, config['MechanicalLossesPanel'], name='MechanicalLossesPanel'),
                     )
        
        for Name, index, panel in zip(['Geometry','State Points','Mass Flow - Valves','Mechanical'],indices,self.panels):
            self.AddPage(panel,Name,imageId=index)
            
        #: A dictionary that maps name to panel 
        self.panels_dict = {panel.Name:panel for panel in self.panels}
    
    def get_config_chunks(self):
        chunks = {}
        for panel in self.panels:
            chunks[panel.name] = panel.get_config_chunk()
        return chunks
        
scroll_yaml=(
"""
family : Scroll Compressor

GeometryPanel:
  Vdisp : 104.8e-6 # Displacement volume / revolution [m^3/rev]
  Vratio : 2.2 # Built-in volume ratio [-]
  ro : 0.005 # Orbiting radius [m]
  t : 0.004 # Scroll wrap thickness [m]
  use_offset : False
  delta_offset : 1e-3 # Offset scroll portion gap width [m]
  phi_fi0 : 0.0 # Inner wrap initial involute angle [rad]
  phi_fis : 3.141 # Inner wrap starting involute angle [rad]
  phi_fos : 0.3 # Outer wrap starting involute angle [rad]
  delta_flank : 15e-6 # Flank gap width [m]
  delta_radial : 15e-6 # Radial gap width [m]
  d_discharge : 0.01 # Discharge port diameter [m]
  inlet_tube_length : 0.3 # Inlet tube length [m]
  inlet_tube_ID : 0.02 # Inlet tube inner diameter [m]
  outlet_tube_length : 0.3 # Outlet tube length [m]
  outlet_tube_ID : 0.02 # Outlet tube inner diameter [m]
  disc_curves:
      type : 2Arc
      r2 : 0

MassFlowPanel:
  sa-s1:
      model : IsentropicNozzle
      options : {Xd : 0.8}
  sa-s2:
      model : IsentropicNozzle
      options : {Xd : 0.8}
  inlet.2-sa:
      model : IsentropicNozzle
      options : {Xd : 0.8}
  d1-dd:
      model : IsentropicNozzle
      options : {Xd : 0.8}
  d2-dd:
      model : IsentropicNozzle
      options : {Xd : 0.8}

MechanicalLossesPanel:
  eta_motor : 0.95 # Motor efficiency [-]
  h_shell : 0.01 # Shell air-side heat transfer coefficient [kW/m^2/K]
  A_shell : 0.040536 # Shell Area [m^2]
  Tamb : 298.0 # Ambient temperature [K]
  mu_oil : 0.0086 # Oil viscosity [Pa-s]
  D_upper_bearing : 0.025 # Upper bearing diameter [m]
  L_upper_bearing : 0.025 # Upper bearing length [m]
  c_upper_bearing : 20e-6 # Upper bearing gap [m]
  D_crank_bearing : 0.025 # Crank bearing diameter [m]
  L_crank_bearing : 0.025 # Crank bearing length [m]
  c_crank_bearing : 20e-6 # Crank bearing gap [m]
  D_lower_bearing : 0.025 # Lower bearing diameter [m]
  L_lower_bearing : 0.025 # Lower bearing length [m]
  c_lower_bearing : 20e-6 # Lower bearing gap [m]
  thrust_friction_coefficient : 0.03 # Thrust bearing friction coefficient [-]
  thrust_ID : 0.08 # Thrust bearing inner diameter [m]
  thrust_OD : 0.3 # Thrust bearing outer diameter [m]
  orbiting_scroll_mass : 2.0 # Orbiting scroll mass [kg]
  L_ratio_bearings : 3.0 # Ratio of lengths for bearings [-]
  HTC : 0.0 # Heat transfer coefficient in scrolls [-]
  journal_tune_factor : 1.0 # Journal loss tune factor [-]
  scroll_plate_thickness : 0.008

StatePanel:
  omega : 377.0 # Rotational speed [rad/s]
  inletState: 
      Fluid : R410A
      T : 283.15 #[K]
      rho : 5.75 #[kg/m^3]
  discharge:
      pratio : 2.0

ParametricPanel:
  structured : True

SolverInputsPanel:
  cycle_integrator: RK45
  integrator_options: {epsRK45 : 1e-7}
  eps_cycle : 0.002 # Cycle-Cycle convergence tolerance (RSSE) [-]
  eps_energy_balance : 0.05 # Energy balance convergence tolerance (RSSE) [-]

Plugin:ScrollInjectionPlugin:
    L = 1

"""
)

def get_defaults():
    return yaml.load(scroll_yaml)