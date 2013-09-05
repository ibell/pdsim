# -*- coding: utf-8 -*-

from __future__ import absolute_import
import wx, yaml, os
from panels import recip_panels, pdsim_panels

family_menu_name = 'Recip Compressor'

# Strings for the script
import_string = 'from PDSim.recip.core import Recip\n'
instantiation_string = 'sim = Recip()\n'
additional_imports_string = 'from PDSim.flow.flow_models import ValveModel\n'

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
        self.panels=(recip_panels.GeometryPanel(self, config['GeometryPanel'],name='GeometryPanel'),
                     pdsim_panels.StateInputsPanel(self, config['StatePanel'], name='StatePanel'),
                     recip_panels.MassFlowPanel(self, config['MassFlowPanel'], name='MassFlowPanel'),
                     recip_panels.MechanicalLossesPanel(self, config['MechanicalLossesPanel'], name='MechanicalLossesPanel'),
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
        
recip_yaml = (
"""
family : Recip Compressor

GeometryPanel:
    piston_diameter : 0.02 # Piston diameter [m]
    piston_length : 0.02 # Piston length [m]
    crank_length : 0.01 # Crank length [m]
    connecting_rod_length : 0.04 # Connecting rod length [m]
    dead_volume_perc : 4.0 #Dead volume percentage [%]
    x_TDC : 0.005 # Distance to piston at TDC [m]
    shell_volume : 100e-6 # Shell volume [\uxb3]

MassFlowPanel:
    d_discharge : 0.0059 # Discharge port diameter [m]
    d_suction : 0.0059 # Suction port diameter [m]
    valve_E : 1.93e+11 # Youngs Modulus [Pa]
    valve_d : 0.007 # Valve diameter [m]
    valve_h : 0.0001532 # Valve thickness [m]
    valve_l : 0.018 # Valve length [m]
    valve_a : 0.014 # Valve distance from anchor [m]
    valve_x_stopper : 0.0018 # Valve distance to stopper [m]
    valve_rho : 8000.0 # Valve metal density [kg/m3]
    valve_C_D : 1.17 # Valve drag coefficient [-]

MechanicalLossesPanel:
    eta_motor : 0.95 # Motor efficiency [-]
    h_shell : 0.01 # Shell air-side heat transfer coefficient [kW/m2/K]
    A_shell : 0.040536 # Shell Area [m2]
    Tamb : 298.0 # Ambient temperature [K]
    mu_oil : 0.0086 # Oil viscosity [Pa-s]
    delta_gap : 2e-05 # Gap width [m]

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
"""
)

def get_defaults():
    return yaml.load(recip_yaml)