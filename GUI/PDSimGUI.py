# -*- coding: latin-1 -*-

import sys,os
sys.path.append(os.path.abspath('panels'))
sys.path.append(os.path.abspath('plugins'))

#Imports from wx package
import wx
from wx.lib.mixins.listctrl import CheckListCtrlMixin, ColumnSorterMixin, ListCtrlAutoWidthMixin
from wx.lib.embeddedimage import PyEmbeddedImage
from wx.lib.wordwrap import wordwrap
wx.SetDefaultPyEncoding('latin-1')

#Provided by python
import codecs
from operator import itemgetter
from math import pi
from Queue import Queue, Empty
from multiprocessing import freeze_support, cpu_count

import time
import textwrap
import cPickle
from ConfigParser import SafeConfigParser
import StringIO
import warnings
import random

#Other packages that are required
import numpy as np
import CoolProp.State as CPState
import yaml
import h5py

#PDSim imports
from PDSim.recip.core import Recip
from PDSim.scroll.core import Scroll
from PDSimLoader import RecipBuilder, ScrollBuilder
from PDSim.plot.plots import PlotNotebook
import PDSim

#PDSim GUI imports
import processes
import pdsim_plugins
import default_configs 
import panels.pdsim_panels as pdsim_panels
import panels.recip_panels as recip_panels
import panels.scroll_panels as scroll_panels
import datatypes

# The path to the home folder that will hold everything
home = os.getenv('USERPROFILE') or os.getenv('HOME')
pdsim_home_folder = os.path.join(home,'.pdsim-temp')
if not os.path.exists(pdsim_home_folder):
    os.mkdir(pdsim_home_folder)
    
class ConfigurationManager(object):
    def __init__(self):
        # Load the config file for the GUI
        if os.path.exists(os.path.join(pdsim_home_folder,'gui_config.yaml')):
            self.config = yaml.load(open(os.path.join(pdsim_home_folder,'gui_config.yaml'),'r'))
        else:
            self.config = {}
            
    def set(self, k, v):
        self.config[k] = v
        yaml.dump(self.config,open(os.path.join(pdsim_home_folder,'gui_config.yaml'),'w'))
        
    def get(self, k, default):
        if k not in self.config:
            return default
        else:
            return self.config[k]

GUIconfig = ConfigurationManager()
    
class IntegratorChoices(wx.Choicebook):
    def __init__(self, parent, **kwargs):
        wx.Choicebook.__init__(self, parent, id = wx.ID_ANY, **kwargs)
    
        # Build the choicebook items
        self.pageEuler=wx.Panel(self)
        self.AddPage(self.pageEuler,'Simple Euler')
        tt = 'Number of steps to be taken by the Euler solver per cycle'
        self.EulerNlabel, self.EulerN = pdsim_panels.LabeledItem(self.pageEuler,
                                                  label="Number of Steps [-]",
                                                  value='7000',
                                                  tooltip = tt)
        sizer=wx.FlexGridSizer(cols=2,hgap=3,vgap=3)
        sizer.AddMany([self.EulerNlabel, self.EulerN])
        self.pageEuler.SetSizer(sizer)
        
        self.pageHeun=wx.Panel(self)
        self.AddPage(self.pageHeun,'Heun')
        tt ='Number of steps to be taken by the Heun solver per cycle'
        self.HeunNlabel, self.HeunN = pdsim_panels.LabeledItem(self.pageHeun,
                                                  label="Number of Steps [-]",
                                                  value='7000',
                                                  tooltip = tt)
        sizer=wx.FlexGridSizer(cols=2,hgap=3,vgap=3)
        sizer.AddMany([self.HeunNlabel, self.HeunN])
        self.pageHeun.SetSizer(sizer)

        tt = """The maximum allowed absolute error per step of the solver"""
        self.pageRK45=wx.Panel(self)
        self.AddPage(self.pageRK45,'Adaptive Runge-Kutta 4/5')
        self.RK45_eps_label, self.RK45_eps = pdsim_panels.LabeledItem(self.pageRK45,
                                                  label="Maximum allowed error per step [-]",
                                                  value='1e-8',
                                                  tooltip = tt)
        sizer=wx.FlexGridSizer(cols=2,hgap=3,vgap=3)
        sizer.AddMany([self.RK45_eps_label, self.RK45_eps])
        self.pageRK45.SetSizer(sizer)
    
    def get_script_chunk(self):
        """
        Returns
        -------
        IC_type : string
            One of 'Euler', 'Heun', 'RK45'
        kwargs : dictionary
            entries to be passed to solve() function
        """
        
        if self.GetSelection() == 0:
            IC_type = 'Euler'
            kwargs = dict(EulerN = int(self.EulerN.GetValue()))
        elif self.GetSelection() == 1:
            IC_type = 'Heun'
            kwargs = dict(HeunN = int(self.HeunN.GetValue()))
        else:
            IC_type = 'RK45'
            kwargs = dict(RK45_eps = int(self.self.RK45_eps.GetValue()))
            
        return IC_type, kwargs
        
class SolverInputsPanel(pdsim_panels.PDPanel):
    
    desc_map = {'eps_cycle' : ('Cycle-cycle convergence criterion','-',0.003),
                'eps_energy_balance' : ('Energy balance convergence criterion','kW',0.01)
                }
    
    def __init__(self, parent, configdict,**kwargs):
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
    
        self.IC = IntegratorChoices(self)
        
        sizer_for_solver_inputs = wx.FlexGridSizer(cols = 2)
        
        annotated_values = []
        self.keys_for_config = []
        
        keys = ['eps_cycle', 'eps_energy_balance']
        annotated_values = self.get_annotated_values(keys, config = configdict)
        
        # Build the items and return the list of annotated GUI objects
        annotated_GUI_objects = self.construct_items(annotated_values, 
                                                     sizer = sizer_for_solver_inputs
                                                     )
        
        # ---------------------------------------------------------------------
        # Register terms in the GUI database
        self.main.register_GUI_objects(annotated_GUI_objects)

        from multiprocessing import cpu_count
        
        sizer_advanced = wx.BoxSizer(wx.VERTICAL)
        self.OneCycle = wx.CheckBox(self, label = "Just run one cycle - not the full solution")
        self.plot_every_cycle = wx.CheckBox(self, label = "Open the plots after each cycle (warning - very annoying but good for debug)")
        
        
        sizer_Ncore_max = wx.BoxSizer(wx.HORIZONTAL)
        
        self.label_Ncore_max = wx.StaticText(self, label = "Maximum number of computational cores to use", )
        self.Ncore_max = wx.SpinCtrl(self, value = "10", )
        self.Ncore_max.SetRange(1, max(1, cpu_count() - 1)) # Ensure that on single-core machines it can still run one core
        self.Ncore_max.SetValue(configdict.get('Ncore_max',1))
        av = datatypes.AnnotatedValue('Ncore_max', self.Ncore_max.GetValue(), 'Maximum number of cores to be used for computation [-]', '-')
        self.main.register_GUI_objects([datatypes.AnnotatedGUIObject(av,self.Ncore_max)])
        
        sizer_Ncore_max.AddMany([self.label_Ncore_max, self.Ncore_max])
        
        sizer_advanced.AddMany([self.OneCycle,
                                self.plot_every_cycle,
                                sizer_Ncore_max])
        
        # Layout the sizers
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(pdsim_panels.HeaderStaticText(self, 'Cycle Integrator Selection'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(self.IC,0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        sizer.Add(pdsim_panels.HeaderStaticText(self, 'Solver Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_solver_inputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        sizer.Add(pdsim_panels.HeaderStaticText(self, 'Advanced and Debug options'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_advanced, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        
        self.SetSizer(sizer)
        sizer.Layout()
    
    def get_config_chunk(self):
  
        IC_type, kwargs = self.IC.get_script_chunk()
        configdict = dict(eps_cycle = self.main.get_GUI_object_value('eps_cycle'),
                          eps_energy_balance = self.main.get_GUI_object_value('eps_energy_balance'),
                          cycle_integrator = IC_type,
                          integrator_options = kwargs,
                          Ncore_max = self.Ncore_max.GetValue()
                          )
        return configdict
        
    def get_script_chunks(self):

        IC_type, kwargs = self.IC.get_script_chunk()
        print IC_type, kwargs
        eps_cycle = self.main.get_GUI_object('eps_cycle').GetValue()
        eps_energy_balance = self.main.get_GUI_object('eps_energy_balance').GetValue()
        
        return textwrap.dedent(
            """
            t1=time.clock()
            sim.connect_callbacks(step_callback = sim.step_callback, 
                                  endcycle_callback = sim.endcycle_callback,
                                  heat_transfer_callback = sim.heat_transfer_callback,
                                  lumps_energy_balance_callback = sim.lump_energy_balance_callback
                                  )
            sim.precond_solve(key_inlet = 'inlet.1',
                              key_outlet = 'outlet.2',
                              pipe_abort = pipe_abort,
                              solver_method = \"{IC_type:s}\",
                              OneCycle = {OneCycle:s},
                              plot_every_cycle = {plot_every_cycle:s},
                              hmin = 1e-8, # hard-coded
                              eps_energy_balance = {eps_energy_balance:s},
                              eps_cycle = {eps_cycle:s},
                              )
            print 'time taken',time.clock()-t1
            """.format(RK_eps = self.IC.RK45_eps.GetValue(),
                       eps_cycle = str(eps_cycle),
                       eps_energy_balance = str(eps_energy_balance),
                       OneCycle = str(self.OneCycle.GetValue()),
                       plot_every_cycle = str(self.plot_every_cycle.GetValue()),
                       IC_type = str(IC_type)
                       )
                   )
        
class SolverToolBook(wx.Toolbook):
    def __init__(self, parent, configdict, id=-1):
        wx.Toolbook.__init__(self, parent, -1, style=wx.BK_LEFT)
        il = wx.ImageList(32, 32)
        indices=[]
        for imgfile in ['Geometry.png','MassFlow.png']:
            ico_path = os.path.join('ico',imgfile)
            indices.append(il.Add(wx.Image(ico_path,wx.BITMAP_TYPE_PNG).ConvertToBitmap()))
        self.AssignImageList(il)
        
        #Make the panels.  Name should be consistent with configuration file
        self.SolverPanel = SolverInputsPanel(self, configdict['SolverInputsPanel'], name = 'SolverInputsPanel')
        self.ParaPanel = pdsim_panels.ParametricPanel(self, configdict['ParametricPanel'], name='ParametricPanel')
        self.panels = (self.SolverPanel, self.ParaPanel)
        
        for Name,index,panel in zip(['Params','Parametric'],indices,self.panels):
            self.AddPage(panel, Name, imageId=index)
            
    def get_config_chunks(self):
        chunks = {}
        for panel in self.panels:
            if hasattr(panel,'get_config_chunk'):
                chunks[panel.name] = panel.get_config_chunk()
        return chunks
    
    def get_script_chunks(self):
        """
        Pull all the values out of the child panels
        
        The values are written into the script file that will be execfile-d
        """
        chunks = []
        for panel in self.panels:
            chunks.append('#############\n# From '+panel.Name+'\n############\n')
            if hasattr(panel,'get_script_params'):
                chunks.append(panel.get_script_params())
            if hasattr(panel,'get_script_chunks'):
                chunks.append(panel.get_script_chunks())
        return chunks

class WriteOutputsPanel(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self,parent)
        
        file_list = ['Temperature', 'Pressure', 'Volume', 'Density','Mass']
        #Create the box
        self.file_list = wx.CheckListBox(self, -1, choices = file_list)
        #Make them all checked
        self.file_list.SetCheckedStrings(file_list)
        
        btn = wx.Button(self,label='Select directory')
        btn.Bind(wx.EVT_BUTTON,self.OnWrite)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.file_list)
        sizer.Add(btn)
        self.SetSizer(sizer)
        
        self.Simulation = None
    
    def set_data(self,Simulation):
        """
        Set the internal simulation data for saving to file
        """
        self.Simulation=Simulation
        
    def OnWrite(self,event):
        """
        Event handler for selection of output folder for writing of files
        """
        dlg = wx.DirDialog(self, "Choose a directory:",
                          style=wx.DD_DEFAULT_STYLE|wx.DD_DIR_MUST_EXIST,
                           #| wx.DD_CHANGE_DIR
                           defaultPath=os.path.abspath(os.curdir)
                           )

        if dlg.ShowModal() == wx.ID_OK:
            self.WriteToFiles(dlg.GetPath())

        # Only destroy a dialog after you're done with it.
        dlg.Destroy()
    
    def WriteToFiles(self,dir_path):
        """
        Write the selected data to files in the folder given by dir_path
        """
        if self.Simulation is None:
            raise ValueError('Simulation data must be provied to WriteOutputsPanel')
        
        outputlist = self.file_list.GetCheckedStrings()
        #List of files that will be over-written
        OWList = [file+'.csv' for file in outputlist if os.path.exists(os.path.join(dir_path,file+'.csv'))]

        if OWList: #if there are any files that might get over-written
            
            dlg = wx.MessageDialog(None,message="The following files will be over-written:\n\n"+'\n'.join(OWList),caption="Confirm Overwrite",style=wx.OK|wx.CANCEL)
            if not dlg.ShowModal() == wx.ID_OK:
                #Don't do anything and return
                return
            
        for file in outputlist:
            if file == 'Pressure':
                xmat = self.Simulation.t
                ymat = self.Simulation.p
                pre = 'p'
            elif file == 'Temperature':
                xmat = self.Simulation.t
                ymat = self.Simulation.T
                pre = 'T'
            elif file == 'Volume':
                xmat = self.Simulation.t
                ymat = self.Simulation.V
                pre = 'V'
            elif file == 'Density':
                xmat = self.Simulation.t
                ymat = self.Simulation.rho
                pre = 'rho'
            elif file == 'Mass':
                xmat = self.Simulation.t
                ymat = self.Simulation.m
                pre = 'm'
            else:
                raise KeyError
            
            #Format for writing (first column is crank angle, following are data)
            joined = np.vstack([xmat,ymat]).T
            
            data_heads = [pre+'['+key+']' for key in self.Simulation.CVs.keys()]
            headers = 'theta [rad],'+ ','.join(data_heads)
            
            def row2string(array):
                return  ','.join([str(dummy) for dummy in array])
            
            rows = [row2string(joined[i,:]) for i in range(joined.shape[0])]
            s = '\n'.join(rows)
            
            #Actually write to file
            print 'writing data to ',os.path.join(dir_path,file+'.csv')
            fp = open(os.path.join(dir_path, file+'.csv'),'w')
            fp.write(headers+'\n')
            fp.write(s)
            fp.close()
            
        print 'You selected: %s\n' % dir_path

class RunToolBook(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self, parent)
        
        # The action buttons
        self.cmdRunOne = wx.Button(self,-1,'  \nRun!\n  ')
        self.cmdRunOne.Bind(wx.EVT_BUTTON, self.GetTopLevelParent().OnStart)
        self.cmdRunParametric = wx.Button(self,-1,'Go To\nParametric\nTable')
        self.cmdRunParametric.Bind(wx.EVT_BUTTON, self.OnGotoParametric)
        self.cmdAbort = wx.Button(self,-1,'Stop\nAll\nRuns')
        self.cmdAbort.Bind(wx.EVT_BUTTON, self.GetTopLevelParent().OnStop)
        
        # Layout the buttons
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(self.cmdRunOne,0)
        hsizer.Add(self.cmdRunParametric,0)
        hsizer.Add(self.cmdAbort,0)
        
        # Define the containers
        splitter = wx.SplitterWindow(self)
        top_pane = wx.Panel(splitter)
        bottom_pane = wx.Panel(splitter)
        
        # The main running log
        self.main_log_ctrl = wx.TextCtrl(top_pane, style = wx.TE_MULTILINE|wx.TE_READONLY)
        
        #The log controls
        nb = wx.Notebook(bottom_pane)
        self.log_ctrls = []
        
        # Make one log box for each thread possible on this machine 
        for thread in range(cpu_count()-1):
            # Make the textbox
            log_ctrl = wx.TextCtrl(nb, style = wx.TE_MULTILINE|wx.TE_READONLY)
            # Add the page to the notebook
            nb.AddPage(log_ctrl,"Thread #"+str(thread+1))
            #Add the log textbox to the list
            self.log_ctrls.append(log_ctrl)
        
        def WriteLog(event=None):
            FD = wx.FileDialog(None,"Log File Name",defaultDir=os.curdir,
                               style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
            if wx.ID_OK == FD.ShowModal():
                fp = open(FD.GetPath(),'w')
                
                # The main log
                fp.write('Main Log\n---------\n')
                fp.write(self.main_log_ctrl.GetValue())
                fp.write('\n\n')
                
                # Print out the logs from each thread
                for i,log in enumerate(self.log_ctrls):
                    fp.write('Log from thread #'+str(i+1)+'\n-------------------\n')
                    fp.write(log.GetValue())
                    fp.write('\n\n')
                fp.close()
            FD.Destroy()
            
        self.write_log_button = wx.Button(bottom_pane, -1, "Write All Logs to a File")
        self.write_log_button.Bind(wx.EVT_BUTTON,WriteLog)
        
        # Layout the top pane
        top_pane_sizer = wx.BoxSizer(wx.VERTICAL)
        top_pane_sizer.Add(wx.StaticText(top_pane,-1,"Main Output Log:"))
        top_pane_sizer.Add(self.main_log_ctrl,1,wx.EXPAND)
        top_pane.SetSizer(top_pane_sizer)
        
        # Layout the bottom pane
        bottom_pane_sizer = wx.BoxSizer(wx.VERTICAL)
        bottom_pane_sizer.Add(nb, 1, wx.EXPAND)
        bottom_pane_sizer.Add(self.write_log_button, 0)
        bottom_pane.SetSizer(bottom_pane_sizer)
        
        # Configure the splitter
        splitter.SetMinimumPaneSize(100)
        splitter.SplitHorizontally(top_pane, bottom_pane, -100)
        
        # Main layout
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(hsizer) # The buttons
        main_sizer.Add(splitter,1, wx.EXPAND)
        self.SetSizer(main_sizer)
        
        top_pane_sizer.Layout()
        bottom_pane_sizer.Layout()
        
        main_sizer.Layout()
    
    def OnGotoParametric(self, event):
        """
        Go to the parametric table
        """
        self.GetTopLevelParent().MTB.SetSelection(1)
        self.GetTopLevelParent().MTB.SolverTB.SetSelection(1)

class FileOutputDialog(wx.Dialog):
    def __init__(self,Simulations, table_string):
        wx.Dialog.__init__(self,None)
        self.Simulations = Simulations
        self.table_string = table_string
        
        #The root directory selector
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(wx.StaticText(self,label="Output Directory:"))
        self.txtDir = wx.TextCtrl(self,value ='.')
        self.txtDir.SetMinSize((200,-1))
        hsizer.Add(self.txtDir,1,wx.EXPAND)
        self.cmdDirSelect = wx.Button(self,label="Select...")
        self.cmdDirSelect.Bind(wx.EVT_BUTTON,self.OnDirSelect)
        hsizer.Add(self.cmdDirSelect)
        
        #The CSV selections
        file_list = ['Temperature', 'Pressure', 'Volume', 'Density','Mass']
        #Create the box
        self.file_list = wx.CheckListBox(self, choices = file_list)
        #Make them all checked
        self.file_list.SetCheckedStrings(file_list)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(hsizer)
        sizer.AddSpacer(10)
        sizer.Add(wx.StaticText(self,label='CSV files:'))
        sizer.Add(self.file_list)
        
        sizer.AddSpacer(10)
        self.chkPickled = wx.CheckBox(self,label='HDF5 data files (Warning! Can be quite large)')
        self.chkPickled.SetToolTipString('Hint: You can use ViTables (search Google) to open the HDF5 files')
        self.chkPickled.SetValue(True)
        sizer.Add(self.chkPickled)
        
        sizer.AddSpacer(10)
        self.chkTable = wx.CheckBox(self,label='Tabular data')
        self.chkTable.SetValue(True)
        sizer.Add(self.chkTable)
        
        self.cmdWrite = wx.Button(self, label = 'Write!')
        self.cmdWrite.Bind(wx.EVT_BUTTON, self.OnWrite)
        sizer.AddSpacer(10)
        sizer.Add(self.cmdWrite)
        self.SetSizer(sizer)
        
    def OnDirSelect(self,event):
        #
        os.chdir(os.curdir)
        dlg = wx.DirDialog(None, "Choose a directory:",
                           defaultPath = os.path.abspath(os.curdir),
                           style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON
                           )
        if dlg.ShowModal() == wx.ID_OK:
            self.txtDir.SetValue(dlg.GetPath())
        dlg.Destroy()
        
    def OnWrite(self, event):
        """
        
        """
        dir_path = self.txtDir.GetValue()
        if not os.path.exists(dir_path):
            dlg = wx.MessageDialog(None, message = 'Selected output directory does not exist.  Please select a folder then try again')
            dlg.ShowModal()
            dlg.Destroy()
            return
        
        for i, sim in enumerate(self.Simulations):
            if (self.file_list.GetCheckedStrings() or 
                self.chkPickled.GetValue()):
                
                run_path = 'RunNumber{0:04d}'.format(i+1)
                if not os.path.exists(os.path.join(dir_path, run_path)):
                    os.mkdir(os.path.join(dir_path, run_path))
                self.write_csv_files(os.path.join(dir_path, run_path), sim)
            
            if self.chkPickled.GetValue():
                self.write_pickle(os.path.join(dir_path, run_path), sim)
            if self.chkTable.GetValue():
                fp = open(os.path.join(dir_path,'ResultsTable.csv'),'w')
                fp.write(self.table_string)
                fp.close()
        self.Destroy()
    
    def write_pickle(self, dir_path, sim):
        from plugins.HDF5_plugin import HDF5Writer
        hdf5_path = os.path.join(dir_path,'Simulation.h5')
        HDF5 = HDF5Writer()
        HDF5.write_to_file(sim, hdf5_path)
        
    def write_csv_files(self, dir_path, sim):
        """
        Write the selected data to files in the folder given by dir_path
        """
        
        outputlist = self.file_list.GetCheckedStrings()
            
        #List of files that will be over-written
        OWList = [file+'.csv' for file in outputlist if os.path.exists(os.path.join(dir_path, file+'.csv'))]

        if OWList: #if there are any files that might get over-written
            
            dlg = wx.MessageDialog(None, message="The following files will be over-written:\n\n"+'\n'.join(OWList),caption="Confirm Overwrite",style=wx.OK|wx.CANCEL)
            if not dlg.ShowModal() == wx.ID_OK:
                #Don't do anything and return
                return wx.ID_CANCEL

        for file in outputlist:
            if file == 'Pressure':
                xmat = sim.t
                ymat = sim.p
                pre = 'p'
            elif file == 'Temperature':
                xmat = sim.t
                ymat = sim.T
                pre = 'T'
            elif file == 'Volume':
                xmat = sim.t
                ymat = sim.V
                pre = 'V'
            elif file == 'Density':
                xmat = sim.t
                ymat = sim.rho
                pre = 'rho'
            elif file == 'Mass':
                xmat = sim.t
                ymat = sim.m
                pre = 'm'
            else:
                raise KeyError
            
            #Format for writing (first column is crank angle, following are data)
            joined = np.vstack([xmat,ymat]).T
            
            data_heads = [pre+'['+key+']' for key in sim.CVs.keys]
            headers = 'theta [rad],'+ ','.join(data_heads)
            
            def row2string(array):
                return  ','.join([str(dummy) for dummy in array])
            
            rows = [row2string(joined[i,:]) for i in range(joined.shape[0])]
            s = '\n'.join(rows)
            
            #Actually write to file
            print 'writing data to ',os.path.join(dir_path,file+'.csv')
            fp = open(os.path.join(dir_path,file+'.csv'),'w')
            fp.write(headers+'\n')
            fp.write(s)
            fp.close()
    
class OutputDataPanel(pdsim_panels.PDPanel):
    def __init__(self, parent, runs, **kwargs):
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        import wx.lib.agw.flatmenu as FM
        from wx.lib.agw.fmresources import FM_OPT_SHOW_CUSTOMIZE, FM_OPT_SHOW_TOOLBAR, FM_OPT_MINIBAR
        
        self._mb = FM.FlatMenuBar(self, wx.ID_ANY, 10, 5, options = FM_OPT_SHOW_TOOLBAR)
        
        self.LoadButton = pdsim_panels.HackedButton(self._mb, label = 'Load Results')
        self._mb.AddControl(self.LoadButton)
        self.LoadButton.Bind(wx.EVT_BUTTON, self.OnLoadRuns)
        
        self.OutputTree = pdsim_panels.OutputTreePanel(self, runs)
        self.AnnotationTarget = wx.StaticText(self, label='Annnotation:')
        
        #Layout of the panel
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self._mb, 0, wx.EXPAND)
        sizer.Add(self.AnnotationTarget,0,wx.EXPAND)
        sizer.Add(self.OutputTree, 1, wx.EXPAND)
        self.SetSizer(sizer)
        sizer.Layout()
    
    def add_runs(self, runs):
        """
        
        Parameters
        ----------
        runs : list
            A list of paths to HDF5 files to be loaded into the GUI
        """
        self.OutputTree.add_runs(runs)
        
    def remove_run(self, index):
        """
        Remove a run at the index given by ``index``
        """
        self.OutputTree.remove_run(index)
        
    def OnLoadRuns(self, event = None):
        """
        Load a HDF5 of output from PDSimGUI
        """
        
        FD = wx.FileDialog(None,"Load Runs",defaultDir = pdsim_home_folder,
                           wildcard = 'PDSim Runs (*.h5)|*.h5',
                           style=wx.FD_OPEN|wx.FD_MULTIPLE|wx.FD_FILE_MUST_EXIST)
        if wx.ID_OK == FD.ShowModal():
            file_paths = FD.GetPaths()
            for file in file_paths:
                self.OutputTree.add_runs(h5py.File(file,'r'))
                print 'added',file
        FD.Destroy()
                
    def OnWriteFiles(self, event):
        """
        Event that fires when the button is clicked to write a selection of things to files
        """
        table_string = self.ResultsList.AsString()
        dlg = FileOutputDialog(self.results, table_string = table_string)
        dlg.ShowModal()
        dlg.Destroy()
        
    def OnRefresh(self, event):
        self.rebuild()
        
class OutputsToolBook(wx.Toolbook):
    def __init__(self, parent, configdict):
        wx.Toolbook.__init__(self, parent, wx.ID_ANY, style=wx.BK_LEFT)
        il = wx.ImageList(32, 32)
        indices=[]
        for imgfile in ['Geometry.png','MassFlow.png']:
            ico_path = os.path.join('ico',imgfile)
            indices.append(il.Add(wx.Image(ico_path,wx.BITMAP_TYPE_PNG).ConvertToBitmap()))
        self.AssignImageList(il)
        
        runs = []

        self.PlotsPanel = wx.Panel(self)
        self.DataPanel = OutputDataPanel(self, runs = runs, name = 'OutputDataPanel')
        
        #Make a  instance
        self.panels = (self.DataPanel, self.PlotsPanel)
        #self.panels = (wx.Panel(self), wx.Panel(self))#self.PlotsPanel)
        for Name,index,panel in zip(['Data','Plots'],indices,self.panels):
            self.AddPage(panel,Name,imageId=index)
            
        self.PN = None
            
    def plot_outputs(self, sim = None):
        # First call there is no plot notebook in existence
        if self.PN is None:
            self.PN = PlotNotebook(sim, self.PlotsPanel)
            sizer = wx.BoxSizer(wx.VERTICAL)
            sizer.Add(self.PN, 1, wx.EXPAND)
            self.PlotsPanel.SetSizer(sizer)
            sizer.Layout()
            self.PlotsPanel.Fit() ##THIS IS VERY IMPORTANT!!!!!!!!!!! :)
        else:
            self.PN.update(sim)
            
class MainToolBook(wx.Toolbook):
    def __init__(self, parent, configdict):
        wx.Toolbook.__init__(self, parent, -1, style=wx.BK_TOP)
        il = wx.ImageList(32, 32)
        indices=[]
        for imgfile in ['Inputs.png','Solver.png','Solver.png','Outputs.png']:
            ico_path = os.path.join('ico',imgfile)
            indices.append(il.Add(wx.Image(ico_path,wx.BITMAP_TYPE_PNG).ConvertToBitmap()))
        self.AssignImageList(il)
        
        # Get the family module to be used
        config_family = configdict['family']
        family = parent.families_dict[config_family]
        
        self.InputsTB = family.InputsToolBook(self, configdict)
        self.SolverTB = SolverToolBook(self, configdict)
        self.RunTB = RunToolBook(self)
        self.OutputsTB = OutputsToolBook(self, configdict)
        
        self.Parent.family_module = family
        
        self.panels=(self.InputsTB,self.SolverTB,self.RunTB,self.OutputsTB)
        for Name,index,panel in zip(['Inputs','Solver','Run','Output'],indices,self.panels):
            self.AddPage(panel,Name,imageId=index)

class MainFrame(wx.Frame):
    def __init__(self, configfile = None, position = None, size = None):
        wx.Frame.__init__(self, None, title = "PDSim GUI", size = (700, 700),style=wx.DEFAULT_FRAME_STYLE |
                          wx.NO_FULL_REPAINT_ON_RESIZE)
        
        self.GUI_object_library = {}
        
        if configfile is None: #No file name or object passed in
            
            defaultconfig = os.path.join('configs','default.pdc')
            
            #First see if a command line option provided
            if '--config' in sys.argv:
                raise NotImplementedError('Not currently accepting command line --config argument')
                i = sys.argv.index('--config')
                _configfile = sys.argv[i+1]
                if os.path.exists(_configfile):
                    self.config = yaml.load(open(_configfile, 'rb'))
                else:
                    warnings.warn('Sorry but your --config file "'+_configfile+'" is not found, loading the default configuration')
                    self.config = yaml.load(open(configfile, 'rb'))
                
            #Then see if there is a file at configs/default.pdc
            elif os.path.exists(defaultconfig):
                self.config = yaml.load(open(defaultconfig,'rb'))
                
            # Then use the internal default scroll
            else:
                self.config = default_configs.get_defaults('scroll')
                
        else:
            # Use the config dictionary passed in
            self.config = configfile
        
        # Get the simulation type (recip, scroll, ...)
        try:
            self.machinefamily = self.config['family']
        except IndexError:
            raise ValueError('configuration file does not have the top-level key "family"')
            
        # The position and size are needed when the frame is rebuilt, but not otherwise
        if position is None:
            position = (-1,-1)
        if size is None:
            size = (-1,-1)
        
        # Use the builder function to rebuild using the configuration objects
        self.build()
        
        # Set up redirection of input and output to logging wx.TextCtrl
        # Taken literally from http://www.blog.pythonlibrary.org/2009/01/01/wxpython-redirecting-stdout-stderr/
        class RedirectText(object):
            def __init__(self,aWxTextCtrl):
                self.out=aWxTextCtrl
            def write(self, string):
                wx.CallAfter(self.out.AppendText, string)
#            def flush(self):
#                return None
                
        redir=RedirectText(self.MTB.RunTB.main_log_ctrl)
        sys.stdout=redir
        sys.stderr=redir
        
        self.SetPosition(position)
        self.SetSize(size)
        
        self.worker = None
        self.workers = None
        self.WTM = None
        
        #: A thread-safe queue for the processing of the results 
        self.results_list = Queue()
        
        # Bind the idle event handler that will always run and
        # deal with the results
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.OnIdle, self.timer)
        self.timer.Start(1000) #1000 ms between checking the queue
        
        self.Bind(wx.EVT_CLOSE, self.OnClose)
    
    def register_GUI_objects(self, annotated_GUI_objects):
        """
        Register the GUI objects in the top-level database
        
        Parameters
        ----------
        annotated_GUI_objects : list of `AnnotatedGUIObject` instances
        """
        # If a single value, convert into a list
        if not isinstance(annotated_GUI_objects,list):
            annotated_GUI_objects = [annotated_GUI_objects]
            
        for o in annotated_GUI_objects:
            if not isinstance(o, datatypes.AnnotatedGUIObject):
                raise TypeError('You can only register lists of AnnotatedGUIObjects. The bad item is:' + str(o))
            if o.key in self.GUI_object_library:
                raise KeyError('Your key [{k:s}] is already in the parameter library'.format(k = o.key))
            if o.annotation in [oo.annotation for oo in self.GUI_object_library.itervalues()]:
                raise KeyError('Your annotation [{a:s}] is already in the parameter library for key [{k:s}]'.format(a = o.annotation, k = o.key))
        
            self.GUI_object_library[o.key] = o
        
    def get_GUI_object(self, key):
        """
        Return an annotated GUI object
        
        Returns
        -------
        o : :class:`AnnotatedGUIObject` instance
        """
        return self.GUI_object_library[key]
    
    def get_GUI_object_value(self, key):
        """
        Return an annotated GUI object's value, converted to a floating point if possible, otherwise as integer or string
        
        Returns
        -------
        val : int,float,string, bool
            value of the GUI object, converted away from string if possible
        """
        val = self.GUI_object_library[key].GetValue()
        if isinstance(val, bool): return val
        try:
            return int(val) #Convert to integer
        except (ValueError,TypeError):
            try:
                return float(val) #Convert to floating point
            except (TypeError,ValueError):
                return val
    
    def get_GUI_object_value_dict(self):
        """
        Return an dictionary of the values of the annotated items in the GUI
        """
        return {key:self.get_GUI_object_value(key) for key in self.get_GUI_object_dict().keys()}
    
    def set_GUI_object_value(self, key, val):
        """
        Set the value of an annotated GUI object's value, converted to the same 
        type as the value of the item is currently
        
        Parameters
        ----------
        val : varied
            Value of the object, can be a class instance to set more complicated
            data types like :class:`CoolProp.State.State` classes
        """
        
        # Get the value of the current target
        current_val = self.GUI_object_library[key].GetValue()
        
        # Get the type of the input to these units 
        type_val = type(current_val)
        
        try:
            # Convert the type of the input to these units
            converted_val = type_val(val)
            
            # Set the item
            self.GUI_object_library[key].SetValue(converted_val)
        
        except TypeError:
            #Try to set the value without converting the value
            self.GUI_object_library[key].SetValue(val)
            
    def get_GUI_object_dict(self):
        """
        Return the dictionary of all the annotated GUI objects
        """
        return self.GUI_object_library
    
    def get_logctrls(self):
        """
        Return a list of the wx.TextCtrl targets for the logs for each thread
        """
        return self.MTB.RunTB.log_ctrls
        
    def rebuild(self, configfile, family_module):
        """
        Destroy everything in the main frame and recreate 
        the contents based on parsing the config file
        """
        # Create a new instance of the MainFrame class using the 
        # new configuration file name and the current location of
        # the frame
        position = self.GetPosition()
        size = self.GetSize()
        try:
            frame = MainFrame(configfile, position=position, size=size)
            frame.Show()
        except:
            raise
        else:
            # Destroy the current MainFrame
            self.Destroy()
            
        self.family_module = family_module
        
    def script_header(self):
        import CoolProp, PDSim
        time_generation = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
        pdsim_version = PDSim.__version__
        coolprop_version = CoolProp.__version__+ ' git revision: '+str(CoolProp.__gitrevision__)
        
        return textwrap.dedent(
            """
            # This file was automatically generated by PDSimGUI on {time_generation:s}
            # Versions used:
            # PDSim: {pdsim_version:s}
            # CoolProp: {coolprop_version:s}
            # 
            # To run this file, in a console, type
            #
            #     python this_file_name.py
            #
            # where this_file_name.py is the name of this file
            #
            # In python 2.7 make a/b always give the double division even 
            # if both a and b are integers - thus 2/3 now yields 0.666666666... 
            # rather than 0
            from __future__ import division
            
            """.format(coolprop_version = coolprop_version,
                       pdsim_version = pdsim_version,
                       time_generation = time_generation)
            )
        
    def script_default_imports(self, plugin_paths):
        
        s = textwrap.dedent(
            """
            # General python imports
            import time, sys, os
            
            # Plotting and numeric things
            from math import pi
            import numpy as np
            from matplotlib import pyplot as plt
            
            # Imports from PDSim that are needed for all families of machines
            from PDSim.flow.flow import FlowPath
            from PDSim.flow.flow_models import IsentropicNozzleWrapper
            from PDSim.core.containers import ControlVolume
            from PDSim.core.core import Tube
            from PDSim.core.motor import Motor
            from PDSim.plot.plots import debug_plots
            
            # Imports from CoolProp
            from CoolProp import State
            from CoolProp import CoolProp as CP
            """
            )
        
        if plugin_paths:
            s += textwrap.dedent(
            """
            # Add the paths for any additional plugin folders (hard-coded absolute paths)
            sys.path.extend({plugin_paths:s})
            """).format(plugin_paths = str(plugin_paths))
        
        return s
        
    def build(self):
        """ Build the entire GUI """
        
        self.make_menu_bar()
        self.load_families(self.FamiliesMenu)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.MTB = MainToolBook(self, self.config)
        sizer.Add(self.MTB, 1, wx.EXPAND)
        self.SetSizer(sizer)
        sizer.Layout()
        
        self.worker=None
        
        self.load_plugins(self.PluginsMenu, self.config)
    
    def build_simulation_script(self, run_index = 1):
        """
        Build the simulation script file
        
        Parameters
        ----------
        run_index : integer, or a string
        
        """
        self.script = []
        
        def indent_chunk(chunks, N):
            if isinstance(chunks,(list,tuple)):
                lines = '\n'.join(chunks)
                lines = lines.split('\n')
            else:
                lines = [chunks]
            
            new_lines = []
            for line in lines:
                new_lines.append('    '*N+line)
            return '\n'.join(new_lines)+'\n'
        
        plugin_chunks = dict(pre_import = '',
                             post_import = '',
                             pre_build = '',
                             post_build = '',
                             pre_build_instantiation = '',
                             post_build_instantiation = '',
                             pre_run = '',
                             post_run = '',
                             )
        
        #Apply any plugins in use - this is the last step in the building process
        if hasattr(self,'plugins_list'):
            for plugin in self.plugins_list:
                if plugin.is_activated():
                    # Get the chunks for this plugin
                    chunks = plugin.get_script_chunks()
                    
                    # Check no non-allowed chunks were returned
                    plugin._check_plugin_chunks(chunks)
                    
                    # Append the keys
                    for key in plugin_chunks.keys():
                        if key in chunks:
                            plugin_chunks[key] += chunks[key]
            
        if hasattr(self,'plugins_list') and any([plugin.is_activated() for plugin in self.plugins_list]):
            #  Get the paths for all the plugin folders if at least one plugin is enabled
            plugin_paths = str(GUIconfig.get('plugin_dirs', default = []))
        else:
            plugin_paths = []
            
        #Get the header for the script
        self.script.append(self.script_header())
        if plugin_chunks['pre_import']:
            self.script.append('####### BEGIN PLUGIN INJECTED CODE (pre-import) ############### \n'+plugin_chunks['pre_import']+'################ END PLUGIN INJECTED CODE ############### \n')
        self.script.append(self.script_default_imports(plugin_paths))
        if hasattr(self,'family_module'):
            self.script.extend(self.family_module.additional_imports_string)
            
        if plugin_chunks['post_import']:
            self.script.append('####### BEGIN PLUGIN INJECTED CODE (post-import) ############### \n'+plugin_chunks['post_import']+'################ END PLUGIN INJECTED CODE ############### \n')
        
        self.script.extend(['def build():\n'])
        if plugin_chunks['pre_build']:
            self.script.extend('####### BEGIN PLUGIN INJECTED CODE (pre-build) ############### \n'+indent_chunk([plugin_chunks['pre_build']],1)+'################ END PLUGIN INJECTED CODE ############### \n')
        self.script.extend([indent_chunk(self.family_module.import_string,1)])
        if plugin_chunks['pre_build_instantiation']:
            self.script.extend('####### BEGIN PLUGIN INJECTED CODE (pre-build-instantiation) ############### \n'+indent_chunk([plugin_chunks['pre_build_instantiation']],1)+'################ END PLUGIN INJECTED CODE ############### \n')
        self.script.extend([indent_chunk(self.family_module.instantiation_string,1)])
        self.script.extend(['    sim.run_index = {run_index:s}\n'.format(run_index = str(run_index))])
        if plugin_chunks['post_build_instantiation']:
            self.script.extend('####### BEGIN PLUGIN INJECTED CODE (post-build-instantiation) ############### \n'+indent_chunk([plugin_chunks['post_build_instantiation']],1)+'################ END PLUGIN INJECTED CODE ############### \n')
        inputs_chunks = self.MTB.InputsTB.get_script_chunks()
        self.script.extend(indent_chunk(inputs_chunks,1))
        if plugin_chunks['post_build']:
            self.script.extend('####### BEGIN PLUGIN INJECTED CODE (post-build) ############### \n'+indent_chunk([plugin_chunks['post_build']],1)+'################ END PLUGIN INJECTED CODE ############### \n')
        self.script.extend(['    return sim\n\n'])

        self.script.extend(['def run(sim, pipe_abort = None):\n'])
        self.script.extend(indent_chunk([plugin_chunks['pre_run']],1))
        solver_chunks = self.MTB.SolverTB.get_script_chunks()
        self.script.extend(indent_chunk(solver_chunks,1))
        self.script.extend(indent_chunk([plugin_chunks['post_run']],1))
        
        self.script.extend(["if __name__ == '__main__':\n    sim = build()\n    run(sim)"])
        
        # Create a string in the form xxxxxxxx where each x is in '0123456789abcdef'
        # This is a simple way to ensure that the file name is unique.  
        random_string = ''.join([random.choice('0123456789abcdef') for i in range(12)])
        
        # Make the full file name
        fName = 'script_' + random_string + '.py'
        
        # Add the PDSim temp folder to the python path
        if pdsim_home_folder not in sys.path:
            sys.path.append(pdsim_home_folder) 
        
        f = open(os.path.join(pdsim_home_folder,fName),'w')
        for chunk in self.script:
            f.write(chunk)
        f.close()

        return fName
    
    def run_batch(self, sims):
        """
        Run a list of simulations
        """
        if self.WTM is None:
            self.MTB.SetSelection(2)
            self.WTM = processes.WorkerThreadManager(sims, 
                                                     self.get_logctrls(),
                                                     done_callback = self.deliver_result,
                                                     main_stdout = self.MTB.RunTB.main_log_ctrl,
                                                     Ncores = self.MTB.SolverTB.SolverPanel.Ncore_max.GetValue()
                                                     )
            self.WTM.setDaemon(True)
            self.WTM.start()
        else:
            dlg = wx.MessageDialog(None,"Batch has already started.  Wait until completion or kill the batch","")
            dlg.ShowModal()
            dlg.Destroy()
            
    def deliver_result(self, hdf5_path = None):
        if hdf5_path is not None:
            self.results_list.put(hdf5_path)
            wx.CallAfter(self.MTB.RunTB.main_log_ctrl.WriteText,'Result queued\n')
     
    def load_plugins(self, PluginsMenu, config):
        """
        Load any plugins into the GUI that are found in plugins folder
        
        It is recommended that all classes and GUI elements relevant to the 
        plugin be included in the given python file
        """
        import glob
        self.plugins_list = []
        
        #  Collect all the .py files in the plugins folder (for standard plugins)
        #  as well as the user-specified plugin folders 
        py_files = []
        for file in glob.glob(os.path.join('plugins','*.py')):
            py_files.append(file)
        
        for directory in GUIconfig.get('plugin_dirs', default = []):
            for file in glob.glob(os.path.join(directory,'*.py')):
                py_files.append(file)
            
        # Hold an old copy of sys.path
        import copy
        old_sys_path = copy.copy(sys.path)
        #  Look at each .py file
        for py_file in py_files:
            #  Get the root filename (/path/to/AAA.py --> AAA)
            path, fname = os.path.split(py_file)
            root, ext = os.path.splitext(fname)
            
            # Hack the path to include the directory
            sys.path = [path]
            
            try:
                #  Try to import the file as a module
                mod = __import__(root)
            except IOError:
                print 'could not import module', root
                continue
            
            for term in dir(mod):
                thing = getattr(mod,term)
                try:
                    #  If it is a plugin class
                    if issubclass(thing, pdsim_plugins.PDSimPlugin):
                        
                        #  Instantiate the plugin
                        plugin = thing()
                        
                        #  Give the plugin a link to the main wx.Frame
                        plugin.set_GUI(self)
                        
                        #  Check if it should be enabled, if not, go to the next plugin
                        if not plugin.should_enable():
                            del plugin
                            continue
                                                
                        #  Append an instance of the plugin to the list of plugins
                        self.plugins_list.append(plugin)
                        
                        #  Create a menu item for the plugin
                        menuItem = wx.MenuItem(PluginsMenu, -1, thing.short_description, "", wx.ITEM_CHECK)
                        PluginsMenu.AppendItem(menuItem)
                        #  Bind the event to activate the plugin
                        self.Bind(wx.EVT_MENU, plugin.activate, menuItem)
                                                
                        #  Check if this type of plugin is included in the config
                        #  file
                        if 'Plugin:'+term in config:
                            #  If it is, activate it and check the element
                            #  in the menu
                            plugin.activate(config = config['Plugin:'+term])
                            menuItem.Check(True)
                        
                except TypeError:
                    pass
        
        # Restore system python path
        sys.path = old_sys_path
   
    def load_families(self, FamiliesMenu):
        """
        Load any machine families into the GUI that are found in families folder
        """
        import glob, pkgutil
        self.plugins_list = []
        
        def load_all_modules_from_dir(dirname, keep_loaded = True):
            mods = []
            for importer, package_name, _ in pkgutil.iter_modules([dirname]):
                full_package_name = '%s.%s' % (dirname, package_name)
                if full_package_name not in sys.modules or keep_loaded:
                    module = importer.find_module(package_name
                                ).load_module(full_package_name)
                    mods.append(module)
            return mods
        
        for mod in load_all_modules_from_dir('families'):
            
            # Get the name of the family
            if hasattr(mod,'family_menu_name'):
                family_menu_name = mod.family_menu_name 
            else:
                raise AttributeError('Family must provide the attribute family_menu_name with the name of the menu item')
            
            # Create a menu item for the plugin
            menuItem = wx.MenuItem(FamiliesMenu, -1, family_menu_name, "", wx.ITEM_CHECK)
            
            # Add the menu item
            FamiliesMenu.AppendItem(menuItem)
            
            # Attach a pointer to the module for this family
            if not hasattr(self,'families_dict'): self.families_dict = {}
            
            # Add both the menu label and module names to the dictionary 
            fname = mod.__name__.rsplit('.',1)[1] #families.scroll --> scroll
            self.families_dict[fname] = mod
            self.families_dict[family_menu_name] = mod
            
            # Bind the event to activate the plugin
            self.Bind(wx.EVT_MENU, self.OnChangeSimFamily, menuItem)
                    
    def make_menu_bar(self):
        #################################
        ####       Menu Bar         #####
        #################################
        
        # Menu Bar
        self.MenuBar = wx.MenuBar()
        
        self.File = wx.Menu()
        self.menuFileOpen = wx.MenuItem(self.File, -1, "Open Config from file...\tCtrl+O", "", wx.ITEM_NORMAL)
        self.menuFileSave = wx.MenuItem(self.File, -1, "Save config to file...\tCtrl+S", "", wx.ITEM_NORMAL)
        self.menuFileFlush = wx.MenuItem(self.File, -1, "Flush out temporary files...", "", wx.ITEM_NORMAL)
        self.menuFileConsole = wx.MenuItem(self.File, -1, "Open a python console", "", wx.ITEM_NORMAL)
        self.menuFileQuit = wx.MenuItem(self.File, -1, "Quit\tCtrl+Q", "", wx.ITEM_NORMAL)
        
        self.File.AppendItem(self.menuFileOpen)
        self.File.AppendItem(self.menuFileSave)
        self.File.AppendItem(self.menuFileFlush)
        self.File.AppendItem(self.menuFileConsole)
        self.File.AppendItem(self.menuFileQuit)
        
        self.MenuBar.Append(self.File, "File")
        self.Bind(wx.EVT_MENU,self.OnOpenConsole,self.menuFileConsole)
        self.Bind(wx.EVT_MENU,self.OnConfigOpen,self.menuFileOpen)
        self.Bind(wx.EVT_MENU,self.OnConfigSave,self.menuFileSave)
        self.Bind(wx.EVT_MENU,self.OnFlushTemporaryFolder,self.menuFileFlush)
        self.Bind(wx.EVT_MENU,self.OnQuit,self.menuFileQuit)
        
        self.FamiliesMenu = wx.Menu()
        self.MenuBar.Append(self.FamiliesMenu, "Families")
        
        self.PluginsMenu = wx.Menu()
        #self.load_plugins(self.PluginsMenu)
        self.MenuBar.Append(self.PluginsMenu, "Plugins")
        self.menuPluginsManage = wx.MenuItem(self.File, -1, "Manage plugin folders...", "", wx.ITEM_NORMAL)
        self.PluginsMenu.AppendItem(self.menuPluginsManage)
        self.Bind(wx.EVT_MENU,self.OnManagePluginFolders,self.menuPluginsManage)
        self.PluginsMenu.AppendSeparator()
        
        self.Solve = wx.Menu()
        self.SolveSolve = wx.MenuItem(self.Solve, -1, "Solve\tF5", "", wx.ITEM_NORMAL)
        self.Solve.AppendItem(self.SolveSolve)
        self.MenuBar.Append(self.Solve, "Solve")
        self.Bind(wx.EVT_MENU, self.OnStart, self.SolveSolve)
        
        self.Help = wx.Menu()
        #self.HelpHelp = wx.MenuItem(self.Help, -1, "Help...\tCtrl+H", "", wx.ITEM_NORMAL)
        self.HelpAbout = wx.MenuItem(self.Help, -1, "About", "", wx.ITEM_NORMAL)
        self.HelpScreenShot = wx.MenuItem(self.Help, -1, "Take screenshot", "", wx.ITEM_NORMAL)        
        self.Help.AppendItem(self.HelpAbout)
        self.Help.AppendItem(self.HelpScreenShot)
        self.MenuBar.Append(self.Help, "Help")
        self.Bind(wx.EVT_MENU, lambda event: self.OnTakeScreenShot(event = None), self.HelpScreenShot)
        self.Bind(wx.EVT_MENU, self.OnAbout, self.HelpAbout)
        
        #Actually set it
        self.SetMenuBar(self.MenuBar)        
        
    ################################
    #         Event handlers       #
    ################################
    
    def OnTakeScreenShot(self, event):
        """ Takes a screenshot of the screen at give pos & size (rect). """
        rect = self.GetRect()
        # see http://aspn.activestate.com/ASPN/Mail/Message/wxpython-users/3575899
        # created by Andrea Gavana
        time.sleep(1)        
        print 'taking screenshot'
        # adjust widths for Linux (figured out by John Torres 
        # http://article.gmane.org/gmane.comp.python.wxpython/67327)
        if sys.platform == 'linux2':
            client_x, client_y = self.ClientToScreen((0, 0))
            border_width = client_x - rect.x
            title_bar_height = client_y - rect.y
            rect.width += (border_width * 2)
            rect.height += title_bar_height + border_width
 
        #Create a DC for the whole screen area
        dcScreen = wx.ScreenDC()
 
        #Create a Bitmap that will hold the screenshot image later on
        #Note that the Bitmap must have a size big enough to hold the screenshot
        #-1 means using the current default colour depth
        bmp = wx.EmptyBitmap(rect.width, rect.height)
 
        #Create a memory DC that will be used for actually taking the screenshot
        memDC = wx.MemoryDC()
 
        #Tell the memory DC to use our Bitmap
        #all drawing action on the memory DC will go to the Bitmap now
        memDC.SelectObject(bmp)
 
        #Blit (in this case copy) the actual screen on the memory DC
        #and thus the Bitmap
        memDC.Blit( 0, #Copy to this X coordinate
                    0, #Copy to this Y coordinate
                    rect.width, #Copy this width
                    rect.height, #Copy this height
                    dcScreen, #From where do we copy?
                    rect.x, #What's the X offset in the original DC?
                    rect.y  #What's the Y offset in the original DC?
                    )
 
        #Select the Bitmap out of the memory DC by selecting a new
        #uninitialized Bitmap
        memDC.SelectObject(wx.NullBitmap)
 
        img = bmp.ConvertToImage()
        fileName = "screenshot.png"
        img.SaveFile(fileName, wx.BITMAP_TYPE_PNG)
        print '...saved as screenshot..png'
        
    def OnManagePluginFolders(self, event):
        plugin_dirs = GUIconfig.get('plugin_dirs', default = [])
        
        class DLG(wx.Dialog):
            def __init__(self,parent,plugin_dirs = [],*args,**kwargs):
                wx.Dialog.__init__(self,parent,size = (510, 400),*args,**kwargs)
        
                #  Build the dialog
                #  Build the dialog
                #  Build the dialog
                
                sizer = wx.BoxSizer(wx.VERTICAL)
                
                sizer.Add(wx.StaticText(self, label='Directories to be searched for plugins'))
                
                self.List = wx.ListBox(self)
                self.List.AppendItems(plugin_dirs)
                self.List.SetMinSize((500,100))
                sizer.Add(self.List)
                
                rsizer = wx.BoxSizer(wx.HORIZONTAL)
                
                self.path_selected = wx.StaticText(self, label = ' Click the select button --------->')
                self.path_selected.SetMinSize((300,-1))
                rsizer.Add(self.path_selected, 0, wx.EXPAND)
                
                self.select_button = wx.Button(self,label='Select...')
                rsizer.Add(self.select_button, 0, wx.EXPAND)
                self.select_button.Bind(wx.EVT_BUTTON,self.OnSelect)
                
                self.add_button = wx.Button(self,label='Add')
                rsizer.Add(self.add_button, 0, wx.EXPAND)
                self.add_button.Bind(wx.EVT_BUTTON,self.OnAdd)
                
                self.remove_button = wx.Button(self,label='Remove Sel.')
                rsizer.Add(self.remove_button, 0, wx.EXPAND)
                self.remove_button.Bind(wx.EVT_BUTTON,self.OnRemoveSelected)
                
                self.ok_button = wx.Button(self,label='Accept')
                rsizer.Add(self.ok_button, 0, wx.EXPAND)
                self.ok_button.Bind(wx.EVT_BUTTON,self.OnAccept)
                
                sizer.Add(rsizer,1)
                
                self.SetSizer(sizer)
                sizer.Layout()
                self.Fit()
                
            def OnAccept(self, event):
                dlg = wx.MessageDialog(None,'Folders to be searched for plugins will be applied at the next startup')
                dlg.ShowModal()
                dlg.Destroy()
                
                self.EndModal(wx.ID_OK)
            
            def OnSelect(self, event):
                dlg = wx.DirDialog(self, "Choose a directory:",
                          style=wx.DD_DEFAULT_STYLE|wx.DD_DIR_MUST_EXIST,
                           #| wx.DD_CHANGE_DIR
                           defaultPath = os.path.abspath(os.curdir)
                           )

                if dlg.ShowModal() == wx.ID_OK:
                    self.path_selected.SetLabel(dlg.GetPath())
                
            def OnAdd(self, event):
                pth = self.path_selected.GetLabel()
                if os.path.exists(pth) and pth not in self.List.GetStrings():
                    self.List.Append(pth)
            
            def OnRemoveSelected(self, event):
                sels = self.List.GetSelections()
                for i in reversed(sorted(sels)):
                    self.List.Delete(i)
        
        dlg = DLG(None, plugin_dirs = plugin_dirs)
        
        if dlg.ShowModal() == wx.ID_OK:
            plugin_dirs = dlg.List.GetStrings()
        
            GUIconfig.set('plugin_dirs', plugin_dirs)
        else:
            print 'closed'
            
        dlg.Destroy()
        
    def OnClose(self, event):
        if hasattr(self,'WTM') and self.WTM is not None and self.WTM.isAlive():
            dlg = wx.MessageDialog(None,"Simulations are running - can't quit")
            dlg.ShowModal()
            dlg.Destroy()
        else:
            event.Skip()
    
    def OnOpenConsole(self, event):
        frm = wx.Frame(None, size = (600,400))
        from wx.py.crust import Crust
        console = Crust(frm, intro = 'Welcome to the debug console within PDSim', locals = locals())
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(console,1,wx.EXPAND)
        frm.SetSizer(sizer)
        frm.Show()
        
    def OnConfigOpen(self,event):
        
        default_dir = GUIconfig.get('config_path', default = 'configs')
        
        FD = wx.FileDialog(None,
                           "Load Configuration file",
                           defaultDir=default_dir,
                           style=wx.FD_OPEN)
        
        if wx.ID_OK==FD.ShowModal():
            
            file_path = FD.GetPath()
            configdict = yaml.load(open(file_path,'r'))
            
            #Now rebuild the GUI using the desired configuration file
            self.rebuild(configdict)
            
            # Write current location to config
            current_path,fname = os.path.split(file_path)
            GUIconfig.set('config_path',current_path)
            
        FD.Destroy()
        
    def OnConfigSave(self,event):
        """
        Write the configuration file
        """
        
        default_dir = GUIconfig.get('config_path', default = 'configs')
        
        FD = wx.FileDialog(None,
                           "Save Configuration file",
                           defaultDir=default_dir,
                           style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        
        if wx.ID_OK == FD.ShowModal():
            #Get the file path
            file_path=FD.GetPath()
            
            print 'Writing configuration file to ', file_path   
            
            #Build the config file entry
            string_list = []
            
            time_generation = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
            #Header information
            header_string_template = textwrap.dedent(
                 """
                 # Generated by PDSimGUI on {time_generation:s}
                 family : {family:s}
                 """
                 )
            terms = dict(time_generation = time_generation,
                         family = self.machinefamily)
            header_string = header_string_template.format(**terms)

            string_list.append(unicode(header_string,'utf-8'))
            
            #Do all the "conventional" panels
            for TB in self.MTB.Children:
                
                #Skip anything that isnt a toolbook
                if not isinstance(TB, wx.Toolbook):
                    continue
                
                if hasattr(TB, 'get_config_chunks'):
                    string_list.append(yaml.dump(TB.get_config_chunks()))
            
            for plugin in self.plugins_list:
                if plugin.is_activated():
                    if hasattr(plugin, 'get_config_chunk'):
                        string_list.append(yaml.dump(plugin.get_config_chunk()))
                    
            fp = codecs.open(file_path,'w',encoding = 'latin-1')
            fp.write(u'\n'.join(string_list))
            fp.close()
            
            # Recreate this frame using the config file to make sure it is 100% the same
            configdict = yaml.load(open(file_path,'r'))
            check_frame = MainFrame(configdict)
            check_dict = check_frame.get_GUI_object_dict()
            my_dict = self.get_GUI_object_dict()
            
            if not sorted(check_dict.keys()) == sorted(my_dict.keys()):
                print 'Not all the keys are the same'
            else:
                for k in check_dict.iterkeys():
                    print check_dict[k].GetValue(),my_dict[k].GetValue()
                    
            # Write current location to config
            current_path,fname = os.path.split(file_path)
            GUIconfig.set('config_path',current_path)
            
        FD.Destroy()
        
    def OnStart(self, event):
        """
        Runs the primary inputs without applying the parametric table inputs
        """
        self.MTB.SetSelection(2)
        script_name = self.build_simulation_script()
        self.run_batch([script_name])
            
    def OnStop(self, event):
        """
        Stop Computation.
        
        Send the abort request to the WorkerThreadManager
        """
        if self.WTM is not None:
            self.WTM.abort()
        
    def OnQuit(self, event):
        self.Destroy()
        wx.Exit()
        
    def OnIdle(self, event):
        """
        Do the things that are needed when the GUI goes idle
        
        This is only run every once in a while (see __init__) for performance-sake 
        """
        
        #Add results from the pipe to the GUI
        if not self.results_list.empty():
            print 'readying to get simulation; ',
            sim = self.results_list.get()
            print 'got a simulation'
            
#            self.MTB.OutputsTB.plot_outputs(sim)
            self.MTB.OutputsTB.DataPanel.add_runs([sim])
#            self.MTB.OutputsTB.DataPanel.rebuild()
            
            #Check whether there are no more results to be processed and threads list is empty
            #This means the manager has completed its work - reset it
            if self.results_list.empty() and not self.WTM.threadsList:
                self.WTM = None
        
        if self.results_list.empty() and self.WTM is not None and not self.WTM.threadsList:
            self.WTM = None
    
    def OnAbout(self, event = None):
        if "unicode" in wx.PlatformInfo:
            wx_unicode = '\nwx Unicode support: True\n'
        else:
            wx_unicode = '\nwx Unicode support: False\n'
        import CoolProp
        info = wx.AboutDialogInfo()
        info.Name = "PDSim GUI"
        info.Version = PDSim.__version__
        info.Copyright = "(C) 2012 Ian Bell"
        info.Description = wordwrap(
            "A graphical user interface for the PDSim model\n\n"+
            "wx version: "+wx.__version__+
            wx_unicode+
            "CoolProp version: "+CoolProp.__version__,
            350, wx.ClientDC(self))
        info.WebSite = ("http://pdsim.sourceforge.net", "PDSim home page")
        info.Developers = [ "Ian Bell", "Craig Bradshaw"]

        # Then we call wx.AboutBox giving it that info object
        wx.AboutBox(info)
        
    def OnChangeSimFamily(self, event):
        """
        Change the machine family of the GUI
        """
        
        for menuItem in self.FamiliesMenu.GetMenuItems():
            # Determine which menu item was clicked
            if event.Id == menuItem.Id: break
        
        # Look up the module based on the name on the GUI menu item
        module = self.families_dict[menuItem.ItemLabel]
        
        # Get its defaut config from the family file and rebuild the GUI using 
        # the default values for this family
        self.rebuild(module.get_defaults(), module)
        
    def OnFlushTemporaryFolder(self, events):
        """
        Event that fires on menu item to flush out temporary files.
        
        Checks to see if temp folder exists, if so, removes it
        """
        import shutil, glob
        
        if os.path.exists(pdsim_home_folder):
            N = len(glob.glob(os.path.join(pdsim_home_folder,'*.*')))
            dlg = wx.MessageDialog(None,'There are '+str(N)+' files in the temporary folder.\n\nPress Ok to remove all the temporary files',style = wx.OK|wx.CANCEL)
            if dlg.ShowModal() == wx.ID_OK:    
                shutil.rmtree(pdsim_home_folder)
                print 'removed the folder',pdsim_home_folder 
            dlg.Destroy()
        else:
            dlg = wx.MessageDialog(None,'Temporary folder does not exist', style = wx.OK)
            dlg.ShowModal()
            dlg.Destroy()

class MySplashScreen(wx.SplashScreen):
    """
    Create a splash screen widget.
    """
    def __init__(self, parent=None):
        # This is a recipe to a the screen.
        # Modify the following variables as necessary.
        img = wx.Image(name = os.path.join("imgs","PDSim_logo.png"))
        width, height = img.GetWidth(), img.GetHeight()
        width *= 0.5
        height *= 0.5
        aBitmap = img.Rescale(width,height).ConvertToBitmap()
        splashStyle = wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT
        splashDuration = 2000 # milliseconds
        # Call the constructor with the above arguments in exactly the
        # following order.
        wx.SplashScreen.__init__(self, aBitmap, splashStyle,
                                 splashDuration, parent)
        self.Bind(wx.EVT_CLOSE, self.OnExit)

        wx.Yield()

    def OnExit(self, evt):
        self.Hide()
        evt.Skip()  # Make sure the default handler runs too...
                    
if __name__ == '__main__':
    # The following line is required to allow cx_Freeze 
    # to package multiprocessing properly.  Must be the first line 
    # after if __name__ == '__main__':
    freeze_support()
    
    app = wx.App(False)
    
    if '--nosplash' not in sys.argv:
        Splash=MySplashScreen()
        Splash.Show()
        time.sleep(2.0)
    
    frame = MainFrame() 
    frame.Show(True) 

    app.MainLoop()
    
    