# -*- coding: utf-8 -*-

from __future__ import absolute_import
import wx, yaml, os
from panels import scroll_panels, pdsim_panels
import numpy as np

family_menu_name = 'Scroll Compressor'

# Strings for the script
import_string = 'from PDSim.scroll.core import Scroll\n'
instantiation_string = 'sim = Scroll()\n'
additional_imports_string = ''

def write_to_xlsx(workbook, runs):
    """
    This provides code to the GUI to add additional terms that are interesting for the scroll compressor
    
    workbook : xlsxwriter.Workbook opened with constant_memory = True
    runs : opened HDF5 files from the runs
    """
    
    ws = workbook.add_worksheet('Pressure profiles')
    
    # Adjust these to adjust the spacing around each block of data
    row_idx = 3
    col_idx = 1
    
    ## Output the matrices
    ## Data has to be written by row because workbook opened with constant_memory = True for writing speed
    
    # Adjust these to adjust the spacing around each block of data
    row_idx = 4
    col_idx = 1
    
    # 3 is the number of empty columns between runs
    offset = 3 + 3
    
    # Header 
    my_col_idx = col_idx
    for run in runs:
        run_index = run.get('run_index').value
        ws.write(row_idx - 3, my_col_idx - 1, 'Run index #'+str(run_index))
        
        if run.get('description'):
            description = run.get('description')
            ws.write(row_idx - 3, my_col_idx - 1, 'Run index #'+str(run_index)+': '+description)
            
        my_col_idx += offset
    
    # Column headers
    my_col_idx = col_idx
    for run in runs:
        ws.write(row_idx-2, my_col_idx - 1, 'theta')
        ws.write(row_idx-2, my_col_idx , 'p1')
        ws.write(row_idx-2, my_col_idx + 1, 'p2')
                
        my_col_idx += offset
        
    datas = []
    maxlen = 0
    for run in runs:
        # Each are stored as 1D array, convert to 2D column matrix
        theta = np.array(run.get('summary/theta_profile').value, ndmin = 2)
        p1 = np.array(run.get('summary/p1_profile').value, ndmin = 2)
        p2 = np.array(run.get('summary/p2_profile').value, ndmin = 2)
    
        datas.append(np.r_[theta, p1, p2].T)
        if np.prod(theta.shape) > maxlen:
            maxlen = np.prod(theta.shape)

    # Data
    for r in range(maxlen):
        my_col_idx = col_idx
        for data in datas:

            if r >= data.shape[0]:
                my_col_idx += offset
                continue
                
            # Theta            
            ws.write(r+row_idx-1, my_col_idx-1,data[r, 0])
            
            if not np.isnan(data[r,1]):
                ws.write(r+row_idx-1, my_col_idx,data[r, 1])
            if not np.isnan(data[r,2]):
                ws.write(r+row_idx-1, my_col_idx+1,data[r, 2])
                    
            my_col_idx += offset    

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