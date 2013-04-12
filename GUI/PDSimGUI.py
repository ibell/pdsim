# -*- coding: latin-1 -*-

#Imports from wx package
import wx
from wx.lib.mixins.listctrl import CheckListCtrlMixin, ColumnSorterMixin, ListCtrlAutoWidthMixin
from wx.lib.embeddedimage import PyEmbeddedImage
from wx.lib.wordwrap import wordwrap
wx.SetDefaultPyEncoding('latin-1')

#Provided by python
import os, sys
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
    
    def set_sim(self, simulation):
        
        if self.GetSelection() == 0:
            simulation.cycle_integrator_type = 'Euler'
            simulation.EulerN = int(self.EulerN.GetValue())
        elif self.GetSelection() == 1:
            simulation.cycle_integrator_type = 'Heun'
            simulation.HeunN = int(self.HeunN.GetValue())
        else:
            simulation.cycle_integrator_type = 'RK45'
            simulation.RK45_eps = float(self.RK45_eps.GetValue())

    def set_from_string(self, config_string):
        """
        config_string will be something like Cycle,Euler,7000 or Cycle,RK45,1e-8
        """
        #Chop off the Cycle,
        config_string = config_string.split(',',1)[1]
        
        SolverType, config = config_string.split(',',1)
        if SolverType == 'Euler':
            self.SetSelection(0)
            self.EulerN.SetValue(config)
        elif SolverType == 'Heun':
            self.SetSelection(1)
            self.HeunN.SetValue(config)
        elif SolverType == 'RK45':
            self.SetSelection(2)
            self.RK45_eps.SetValue(config)
        
    def save_to_string(self):
        if self.GetSelection() == 0:
            return 'Cycle = Cycle,Euler,'+self.EulerN.GetValue()
        elif self.GetSelection() == 1:
            return 'Cycle = Cycle,Heun,'+self.HeunN.GetValue()
        else:
            return 'Cycle = Cycle,RK45,'+self.RK45_eps.GetValue()
        
class SolverInputsPanel(pdsim_panels.PDPanel):
    
    desc_map = {'eps_cycle' : ('Cycle-cycle convergence criterion','-')}
    
    def __init__(self, parent, configdict,**kwargs):
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
    
        self.IC = IntegratorChoices(self)
        
        sizer_for_solver_inputs = wx.BoxSizer(wx.HORIZONTAL)
        
        annotated_values = []
        # Loop over the first group of inputs
        for key in ['eps_cycle']:
            # Get the annotation and the units for the term 
            annotation, units = self.desc_map[key]
            # Add the annotated object to the list of objects
            annotated_values.append(datatypes.AnnotatedValue(key, configdict[key], annotation, units))
        
        # Build the items and return the list of annotated GUI objects
        annotated_GUI_objects = self.construct_items(annotated_values, 
                                                     sizer = sizer_for_solver_inputs
                                                     )
        
        # ---------------------------------------------------------------------
        # Register terms in the GUI database
        self.main.register_GUI_objects([annotated_GUI_objects])

        sizer_advanced = wx.BoxSizer(wx.VERTICAL)
        self.OneCycle = wx.CheckBox(self, label = "Just run one cycle - not the full solution")
        self.plot_every_cycle = wx.CheckBox(self, label = "Open the plots after each cycle (warning - very annoying but good for debug)")
        sizer_advanced.AddMany([self.OneCycle,
                                self.plot_every_cycle])
        
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
        sizer.Add(pdsim_panels.HeaderStaticText(self, 'Advanced & Debug options'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_advanced, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        
        self.SetSizer(sizer)
        sizer.Layout()
    
    def get_config_chunk(self):
  
        configdict = dict(eps_cycle = self.main.get_GUI_object_value('eps_cycle'),
                          cycle_integrator = 'RK45',
                          integrator_options = dict(RK_eps = float(self.IC.RK45_eps.GetValue()))
                          )
        return configdict
        
    def get_script_chunks(self):
        
        eps_cycle = self.main.get_GUI_object('eps_cycle').GetValue()
        
        return textwrap.dedent(
            """
            sim.RK45_eps = {RK_eps:s}
            t1=time.clock()
            sim.connect_callbacks(step_callback = sim.step_callback, 
                                  endcycle_callback = sim.endcycle_callback,
                                  heat_transfer_callback = sim.heat_transfer_callback,
                                  lumps_energy_balance_callback = sim.lump_energy_balance_callback
                                  )
            sim.precond_solve(key_inlet='inlet.1',
                              key_outlet='outlet.2',
                              pipe_abort = pipe_abort,
                              solver_method = 'RK45',
                              UseNR = False, #Don't use Newton-Raphson ND solver to determine the initial state
                              OneCycle = {OneCycle:s},
                              plot_every_cycle = {plot_every_cycle:s},
                              hmin = 1e-8 # hard-coded
                              )
            print 'time taken',time.clock()-t1
            """.format(RK_eps = self.IC.RK45_eps.GetValue(),
                       eps_cycle = str(eps_cycle),
                       OneCycle = str(self.OneCycle.GetValue()),
                       plot_every_cycle = str(self.plot_every_cycle.GetValue())
                       )
                   )
            
    def post_get_from_configfile(self, key, config_string):
        """
        Build the integrator chooser 
        
        This function will be called by PDPanel.get_from_configfile
        """
        if key == 'Cycle':
            self.IC.set_from_string(config_string)
        
    def post_prep_for_configfile(self):
        return self.IC.save_to_string()+'\n'
    
    def supply_parametric_term(self):
        pass
        
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
        pane1 = SolverInputsPanel(self, configdict['SolverInputsPanel'], name = 'SolverInputsPanel')
        pane2 = pdsim_panels.ParametricPanel(self, configdict['ParametricPanel'], name='ParametricPanel')
        self.panels = (pane1, pane2)
        
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
        Pull all the values out of the child panels, using the values in 
        self.items and the function get_script_chunks if the panel implements
        it
        
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
        self.cmdRunParametric = wx.Button(self,-1,'Run\nParametric\nTable')
        #self.cmdRunParametric.Bind(wx.EVT_BUTTON, self.GetTopLevelParent().OnStart)
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
        
class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent, ID = wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        ListCtrlAutoWidthMixin.__init__(self)
        
class ResultsList(wx.Panel, ColumnSorterMixin):
    def __init__(self, parent, headers, values, results):
        """
        
        parent : wx.Window
            parent of the Panel
            
        headers: a list of strings
            Each element is the string that will be the header of the column
            
        values: a list of list of values.  
            Each entry in the list should be as long as the number of headers
            
        results : PDSimCore instances
            The simulation runs
        """
        wx.Panel.__init__(self, parent)
        
        #: The list of strings of the header
        self.headers = list(headers)
        #: The values in the table
        self.values = list(values)
        #: The PDSimCore instances that have all the data
        self.results = list(results)
        
        self.list = AutoWidthListCtrl(self, 
                                      style=wx.LC_REPORT | wx.BORDER_NONE
                                      )
        #Build the headers
        for i, header in enumerate(headers):
            self.list.InsertColumn(i, header)
        
        #Add the values one row at a time
        self.itemDataMap = {}
        for i, row in enumerate(self.values):
            #Add an entry to the data map
            self.itemDataMap[i] = tuple(row)
            
            self.list.InsertStringItem(i,str(row[0]))
            self.list.SetItemData(i,i)
            
            for j in range(1,len(row)):
                self.list.SetStringItem(i,j,str(row[j]))
        
        total_width = 0    
        for i in range(len(headers)):
            self.list.SetColumnWidth(i, wx.LIST_AUTOSIZE_USEHEADER)
            total_width += self.list.GetColumnWidth(i)
            
        width_available = self.Parent.GetSize()[0]
        self.list.SetMinSize((width_available,200))
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.list,1,wx.EXPAND)
        
        self.il = wx.ImageList(16, 16)
        SmallUpArrow = PyEmbeddedImage(
            "iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABHNCSVQICAgIfAhkiAAAADxJ"
            "REFUOI1jZGRiZqAEMFGke2gY8P/f3/9kGwDTjM8QnAaga8JlCG3CAJdt2MQxDCAUaOjyjKMp"
            "cRAYAABS2CPsss3BWQAAAABJRU5ErkJggg==")
        SmallDnArrow = PyEmbeddedImage(
            "iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABHNCSVQICAgIfAhkiAAAAEhJ"
            "REFUOI1jZGRiZqAEMFGke9QABgYGBgYWdIH///7+J6SJkYmZEacLkCUJacZqAD5DsInTLhDR"
            "bcPlKrwugGnCFy6Mo3mBAQChDgRlP4RC7wAAAABJRU5ErkJggg==")
        
        self.sm_up = self.il.Add(SmallUpArrow.GetBitmap())
        self.sm_dn = self.il.Add(SmallDnArrow.GetBitmap())
        self.list.SetImageList(self.il, wx.IMAGE_LIST_SMALL)

        ColumnSorterMixin.__init__(self,len(headers)+1)
        
        self.SetSizer(sizer)
        self.SetAutoLayout(True)
        
    def OnSortOrderChanged(self, *args, **kwargs):
        """
        Overload the base class method to resort the internal structures 
        when the table is sorted
        """
        self._sort_objects()
        return ColumnSorterMixin.OnSortOrderChanged(self, *args, **kwargs)
        
    def __getitem__(self, index):
        """
        Provided to be able to index the class
        returns the index of the run returned
        """
        
        return self.results[index]
    
    def _sort_objects(self):
        """
        Sort the internal data structures based on the table sort state
        """
        
        # Sort the output csv table in the same way as the listctrl
        iCol, direction = self.GetSortState()
        
        #If sorted, sort the variables
        if iCol >= 0:
        
            # Get a sorted version of self.values sorted by the column used in list
            values_results = zip(self.values, self.results)
            
            # Sort the results and the rows together
            sorted_things = sorted(values_results, key=itemgetter(iCol))
            
            # Unpack
            self.values, self.results = zip(*sorted_things)
            
            # tuples --> list
            self.values = list(self.values)
            self.results = list(self.results)
        
    def remove_item(self, index):
        """
        Remove the item from the ResultsList instance
        """
        #Remove the item from the data map
        self.itemDataMap.pop(index)
        #Remove the item from the values
        del self.values[index]
        #Remove the item from the results
        del self.results[index]
        
    def get_results(self):
        """
        Return the list of PDSimCore instances
        """
        return self.results
 
    def GetListCtrl(self):
        """
        Required method for ColumnSorterMixin
        
        Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py
        """
        return self.list
    
    def GetSortImages(self):
        """
        Required method for ColumnSorterMixin
        
        Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py
        """
        return (self.sm_dn, self.sm_up)
    
    def AsString(self):
        """
        Return a csv formatted table of the ResultsList
        
        """
        
        #Sort the internal data structures based on table sort
        self._sort_objects()
        
        header_string = [','.join(self.headers)]
        def tostr(row):
            return [str(r) for r in row]
        rows_string = [','.join(tostr(row)) for row in self.values]
        return '\n'.join(header_string+rows_string)

class ColumnSelectionDialog(wx.Dialog):
    def __init__(self, parent, col_options, cols_selected):
        wx.Dialog.__init__(self,parent,size = (800,350))
        
        self.col_options = col_options

        self.selected = [col_options[col] for col in cols_selected]
        self.not_selected = [col_options[col] for col in col_options if col not in cols_selected]
        
        self.col_library_label = wx.StaticText(self, label = 'Available columns:')
        self.col_used_label = wx.StaticText(self, label = 'Selected columns:')
        self.col_library = wx.ListBox(self, choices = self.not_selected, style = wx.LB_EXTENDED)
        self.col_used = wx.ListBox(self, choices = self.selected, style = wx.LB_EXTENDED)
        self.col_library.SetMinSize((300,300))
        self.col_used.SetMinSize((300,300))
        
        #The central column with add and remove buttons
        self.AddAllButton=wx.Button(self, label='All ->')
        self.RemoveAllButton=wx.Button(self, label='<- All')
        self.AddButton=wx.Button(self, label='-->')
        self.RemoveButton=wx.Button(self, label='<--')
        self.AddButton.Bind(wx.EVT_BUTTON,self.OnAdd)
        self.RemoveButton.Bind(wx.EVT_BUTTON,self.OnRemove)
        self.AddAllButton.Bind(wx.EVT_BUTTON,self.OnAddAll)
        self.RemoveAllButton.Bind(wx.EVT_BUTTON,self.OnRemoveAll)
        vsizer = wx.BoxSizer(wx.VERTICAL)
        vsizer.AddMany([self.AddAllButton, self.RemoveAllButton])
        vsizer.AddSpacer(40)
        vsizer.AddMany([self.AddButton, self.RemoveButton])

        #The far-right column with up,down, ok, cancel buttons      
        self.Up = wx.Button(self, label='Move Up')
        self.Up.Bind(wx.EVT_BUTTON,self.OnUp)
        self.Down = wx.Button(self, label='Move Down')
        self.Down.Bind(wx.EVT_BUTTON,self.OnDown)
        self.OkButton = wx.Button(self, label='Ok')
        self.OkButton.Bind(wx.EVT_BUTTON,self.OnAccept)
        self.CancelButton = wx.Button(self, label='Cancel')
        self.CancelButton.Bind(wx.EVT_BUTTON,self.OnClose)
        vsizer2 = wx.BoxSizer(wx.VERTICAL)
        vsizer2.AddMany([self.Up,self.Down])
        vsizer2.AddSpacer(40)
        vsizer2.AddMany([self.CancelButton, self.OkButton])
        
        #Layout the dialog
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        vsizer0 = wx.BoxSizer(wx.VERTICAL)
        vsizer0.Add(self.col_library_label)
        vsizer0.Add(self.col_library, 1, wx.EXPAND)
        sizer.Add(vsizer0)
        sizer.AddSpacer(10)
        sizer.Add(vsizer,0,wx.ALIGN_CENTER_VERTICAL)
        sizer.AddSpacer(10)
        vsizer20 = wx.BoxSizer(wx.VERTICAL)
        vsizer20.Add(self.col_used_label)
        vsizer20.Add(self.col_used, 1, wx.EXPAND)
        sizer.Add(vsizer20)
        sizer.AddSpacer(10)
        sizer.Add(vsizer2,0,wx.ALIGN_CENTER_VERTICAL)
        self.SetSizer(sizer)
        sizer.Layout()
        
        #Bind a key-press event to all objects to get Esc key press
        children = self.GetChildren()
        for child in children:
            child.Bind(wx.EVT_KEY_UP,  self.OnKeyPress) 

    def OnKeyPress(self,event):
        """ cancel if Escape key is pressed """
        event.Skip()
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self.EndModal(wx.ID_CANCEL)
        elif event.GetKeyCode() == wx.WXK_RETURN:
            self.EndModal(wx.ID_OK)
        
    def label2attr(self,label):
        for col in self.col_options:
            if self.col_options[col] == label:
                return col
        raise KeyError
        
    def OnAccept(self, event):
        self.EndModal(wx.ID_OK)
        
    def OnClose(self,event):
        self.EndModal(wx.ID_CANCEL)
        
    def OnAddAll(self, event):
        self.selected += self.not_selected
        self.not_selected = []
        self.col_library.SetItems(self.not_selected)
        self.col_used.SetItems(self.selected)
        
    def OnRemoveAll(self, event):
        self.not_selected += self.selected
        self.selected = []
        self.col_library.SetItems(self.not_selected)
        self.col_used.SetItems(self.selected)
        
    def OnAdd(self, event):
        indices = self.col_library.GetSelections()
        labels = [self.col_library.GetString(index) for index in indices]

        for label in reversed(labels):
            i = self.not_selected.index(label)
            self.selected.append(self.not_selected.pop(i))
        self.col_library.SetItems(self.not_selected)
        self.col_used.SetItems(self.selected)
        
    def OnRemove(self, event):
        indices = self.col_used.GetSelections()
        labels = [self.col_used.GetString(index) for index in indices]

        for label in reversed(labels):
            i = self.selected.index(label)
            self.not_selected.append(self.selected.pop(i))
        self.col_library.SetItems(self.not_selected)
        self.col_used.SetItems(self.selected)
        
    def OnUp(self, event):
        indices = self.col_used.GetSelections()
        labels = [self.col_used.GetString(index) for index in indices]
        for label in labels:
            i = self.selected.index(label)
            if i>0:
                #swap item and the previous item
                self.selected[i-1],self.selected[i]=self.selected[i],self.selected[i-1]
        self.col_used.SetItems(self.selected)
        if len(labels) == 1:
            self.col_used.SetSelection(indices[0]-1)
    
    def OnDown(self, event):
        indices = self.col_used.GetSelections()
        labels = [self.col_used.GetString(index) for index in indices]
        for label in labels:
            i = self.selected.index(label)
            if i<len(self.selected)-1:
                #swap item and the next item
                self.selected[i+1],self.selected[i]=self.selected[i],self.selected[i+1]
        self.col_used.SetItems(self.selected)
        if len(labels) == 1:
            self.col_used.SetSelection(indices[0]+1)
    
    def GetSelections(self):
        labels = self.col_used.GetStrings()
        attrs = [self.label2attr(label) for label in labels]
        return attrs

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
        
        self.RefreshButton = pdsim_panels.HackedButton(self._mb, label = 'Refresh')
        self._mb.AddControl(self.RefreshButton)
        #self.BuildButton.Bind(wx.EVT_BUTTON, self.OnBuildTable)
        self.RefreshButton.Disable()
        
        self.RunButton = pdsim_panels.HackedButton(self._mb, label = 'Run Table')
        self._mb.AddControl(self.RunButton)
        #self.RunButton.Bind(wx.EVT_BUTTON, self.OnRunTable)
        self.RunButton.Disable()
        
        self.OutputTree = pdsim_panels.OutputTreePanel(self, runs)
        
        #Layout of the panel
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self._mb, 0, wx.EXPAND)
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
    
    def OnPlotSelected(self, event):
        list_ = self.ResultsList.GetListCtrl()
        
        indices = []
        index = list_.GetFirstSelected()
        sim = self.results[index]
        self.Parent.plot_outputs(sim)
        
        if list_.GetNextSelected(index) != -1:
            dlg = wx.MessageDialog(None,'Sorry, only the first selected row will be used')
            dlg.ShowModal()
            dlg.Destroy()
                
    def OnRemoveSelected(self, event):
        list_ = self.ResultsList.GetListCtrl()
        
        indices = []
        index = list_.GetFirstSelected()
        while index != -1:
            indices.append(index)
            index = list_.GetNextSelected(index)
            
        #Some runs to delete
        if indices:
            #Warn before removing
            dlg = wx.MessageDialog(None,'You are about to remove '+str(len(indices))+' runs.  Ok to confirm', style = wx.OK|wx.CANCEL)
            if dlg.ShowModal() == wx.ID_OK:
                for index in reversed(indices):
                    #Remove the item from the ResultsList
                    self.ResultsList.remove_item(index)
                    
                #Update our copy of the results
                self.results = self.ResultsList.get_results()
                
            dlg.Destroy()
            
            #Rebuild the ResultsList
            self.rebuild()
                
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
        
        # Load the runs into memory
        runa = h5py.File('runa.h5')#, driver = 'core', backing_store = False)
        runb = h5py.File('runb.h5')#, driver = 'core', backing_store = False)
        
        runs = [runa,runb]
#        variables = self.Parent.InputsTB.collect_parametric_terms()
#        self.PlotsPanel = wx.Panel(self)
        self.DataPanel = OutputDataPanel(self, runs = runs, name = 'OutputDataPanel')
        
        #Make a  instance
        self.panels = (self.DataPanel, wx.Panel(self))#self.PlotsPanel)
        #self.panels = (wx.Panel(self), wx.Panel(self))#self.PlotsPanel)
        for Name,index,panel in zip(['Data','Plots'],indices,self.panels):
            self.AddPage(panel,Name,imageId=index)
            
        self.PN = None
            
    def plot_outputs(self, recip = None):
        parent = self.PlotsPanel
        # First call there is no plot notebook in existence
        if self.PN is None:
            self.PN = PlotNotebook(recip,parent)
            sizer = wx.BoxSizer(wx.VERTICAL)
            sizer.Add(self.PN,1,wx.EXPAND)
            parent.SetSizer(sizer)
            parent.Fit() ##THIS IS VERY IMPORTANT!!!!!!!!!!! :)
        else:
            self.PN.update(recip)
            
    def add_output_terms(self, items):
        self.DataPanel.add_output_terms(items)
        
    def change_output_terms(self, key_dict):
        self.DataPanel.change_output_terms(key_dict)
            
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
        
        #Get the simulation type (recip, scroll, ...)
        try:
            self.machinefamily = self.config['family']
        except IndexError:
            raise ValueError('configuration file does not have the top-level key "family"')
            
        #The position and size are needed when the frame is rebuilt, but not otherwise
        if position is None:
            position = (-1,-1)
        if size is None:
            size = (-1,-1)
        
        #Use the builder function to rebuild using the configuration objects
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
                
#        redir=RedirectText(self.MTB.RunTB.log_ctrl)
#        sys.stdout=redir
#        sys.stderr=redir
        
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
                raise TypeError('You can only register lists of AnnotatedGUIObjects')
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
    
    def update_parametric_terms(self):
        """
        Actually update the parametric terms in the parametric table options
        """
        para_terms = self.collect_parametric_terms()
        self.MTB.SolverTB.update_parametric_terms(para_terms)
        
    def collect_parametric_terms(self):
        """
        This function is called to find all the parametric terms that are
        required.
        
        They can be recursively found in:
        - self.items in PDPanel instances
        - collect_parametric_terms in PDPanel instances
        - 
        
        """
        terms = []
        #Loop over the toolbooks and allow them to collect their own terms
        for child in self.MTB.Children:
            if isinstance(child,wx.Toolbook) and hasattr(child,'collect_parametric_terms'):
                terms += child.collect_parametric_terms()
        return terms
        
    def rebuild(self, configfile):
        """
        Destroy everything in the main frame and recreate 
        the contents based on parsing the config file
        """
        # Create a new instance of the MainFrame class using the 
        # new configuration file name and the current location of
        # the frame
        position = self.GetPosition()
        size = self.GetSize()
        frame = MainFrame(configfile, position=position, size=size)
        frame.Show()
        
        #Destroy the current MainFrame
        self.Destroy()
        
    def script_header(self):
        import CoolProp, PDSim
        time_generation = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
        pdsim_version = PDSim.__version__
        coolprop_version = CoolProp.__version__+ ' svn revision: '+str(CoolProp.__svnrevision__)
        
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
            from PDSim.plot.plots import debug_plots
            from PDSim.core.motor import Motor
            
            # Imports from CoolProp
            from CoolProp import State
            from CoolProp import CoolProp as CP
            
            
            """.format(coolprop_version = coolprop_version,
                       pdsim_version = pdsim_version,
                       time_generation = time_generation)
            )
        
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
        
#        self.load_plugins(self.PluginsMenu)
#        #After loading plugins, try to set the parameters in the parametric table
#        self.MTB.SolverTB.flush_parametric_terms()
#        self.MTB.SolverTB.set_parametric_terms()
    
    def build_simulation_script(self, run_index = 1):
        """
        Build the simulation script file
        
        Parameters
        ----------
        run_index : integer, or a string
        
        """
        self.script_chunks = []
        
        #Get the header for the script
        self.script_chunks.append(self.script_header())
        
        def indent_chunk(chunks, N):
            lines = '\n'.join(chunks)
            lines = lines.split('\n')
            
            new_lines = []
            for line in lines:
                new_lines.append('    '*N+line)
            return '\n'.join(new_lines)+'\n'
        
#        #Apply any plugins in use - this is the last step in the building process
#        if hasattr(self,'plugins_list'):
#            for plugin in self.plugins_list:
#                if plugin.is_activated():
#                    raise NotImplementedError
#                    plugin.get_script_chunks()
            
        self.script_chunks.extend(['def build():\n'])
        self.script_chunks.extend(['    from PDSim.scroll.core import Scroll\n    sim = Scroll()\n'])
        self.script_chunks.extend(['    sim.run_index = {run_index:s}\n'.format(run_index = str(run_index))])
        inputs_chunks = self.MTB.InputsTB.get_script_chunks()
        self.script_chunks.extend(indent_chunk(inputs_chunks,1))
        self.script_chunks.extend(['    return sim\n\n'])

        self.script_chunks.extend(['def run(sim, pipe_abort = None):\n'])
        solver_chunks = self.MTB.SolverTB.get_script_chunks()
        self.script_chunks.extend(indent_chunk(solver_chunks,1))
        
        self.script_chunks.extend(["if __name__ == '__main__':\n    sim = build()\n    run(sim)"])
        
        # Create a string in the form xxxxxxxx where each x is in '0123456789abcdef'
        # This is a simple way to ensure that the file name is unique.  
        random_string = ''.join([random.choice('0123456789abcdef') for i in range(12)])
        
        # Make the full file name
        fName = 'script_' + random_string + '.py'
        
        # Add the PDSim temp folder to the python path
        if pdsim_home_folder not in sys.path:
            sys.path.append(pdsim_home_folder) 
        
        f = open(os.path.join(pdsim_home_folder,fName),'w')
        for chunk in self.script_chunks:
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
                                                     main_stdout = self.MTB.RunTB.main_log_ctrl)
            self.WTM.setDaemon(True)
            self.WTM.start()
        else:
            dlg = wx.MessageDialog(None,"Batch has already started.  Wait until completion or kill the batch","")
            dlg.ShowModal()
            dlg.Destroy()
            
    def deliver_result(self, sim = None):
        if sim is not None:
            self.results_list.put(sim)
            wx.CallAfter(self.MTB.RunTB.main_log_ctrl.WriteText,'Result queued\n')
     
    def load_plugins(self, PluginsMenu):
        """
        Load any plugins into the GUI that are found in plugins folder
        
        It is recommended that all classes and GUI elements relevant to the 
        plugin be included in the given python file
        """
        import glob
        self.plugins_list = []
        #Look at each .py file in plugins folder
        for py_file in glob.glob(os.path.join('plugins','*.py')):
            #Get the root filename (/path/to/AAA.py --> AAA)
            fname = py_file.split(os.path.sep,1)[1].split('.')[0]
            
            mods = __import__('plugins.'+fname)
            #Try to import the file as a module
            mod = getattr(mods,fname)
            for term in dir(mod):
                thing = getattr(mod,term)
                try:
                    #If it is a plugin class
                    if issubclass(thing, pdsim_plugins.PDSimPlugin):
                        
                        #Instantiate the plugin
                        plugin = thing()
                        
                        #Give the plugin a link to the main wx.Frame
                        plugin.set_GUI(self)
                        
                        #Check if it should be enabled, if not, go to the next plugin
                        if not plugin.should_enable():
                            del plugin
                            continue
                                                
                        #Append an instance of the plugin to the list of plugins
                        self.plugins_list.append(plugin)
                        
                        #Create a menu item for the plugin
                        menuItem = wx.MenuItem(self.Type, -1, thing.short_description, "", wx.ITEM_CHECK)
                        PluginsMenu.AppendItem(menuItem)
                        #Bind the event to activate the plugin
                        self.Bind(wx.EVT_MENU, plugin.activate, menuItem)
                                                
                        # Check if this type of plugin is included in the config
                        # file
                        for section in self.config_parser.sections():
                            if (section.startswith('Plugin')
                                    and section.split(':')[1] ==  term):
                                # If it is, activate it and check the element
                                # in the menu
                                plugin.activate()
                                menuItem.Check(True)
                                
                                # Pass the section along to the plugin
                                items = self.config_parser.items(section)
                                plugin.build_from_configfile_items(items)
                        
                except TypeError:
                    pass
        
        # Update the parametric terms in the parametric tables because the 
        # plugins might have added terms if they are activated from the config
        # file
        self.update_parametric_terms()
   
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
        
#        self.Bind(wx.EVT_MENU,self.OnChangeSimType,self.TypeScroll)
#        self.Bind(wx.EVT_MENU,self.OnChangeSimType,self.TypeRecip)
        
        self.Solve = wx.Menu()
        self.SolveSolve = wx.MenuItem(self.Solve, -1, "Solve\tF5", "", wx.ITEM_NORMAL)
        self.Solve.AppendItem(self.SolveSolve)
        self.MenuBar.Append(self.Solve, "Solve")
        self.Bind(wx.EVT_MENU, self.OnStart, self.SolveSolve)
        
        self.Help = wx.Menu()
        #self.HelpHelp = wx.MenuItem(self.Help, -1, "Help...\tCtrl+H", "", wx.ITEM_NORMAL)
        self.HelpAbout = wx.MenuItem(self.Help, -1, "About", "", wx.ITEM_NORMAL)
        #self.Help.AppendItem(self.HelpHelp)
        self.Help.AppendItem(self.HelpAbout)
        self.MenuBar.Append(self.Help, "Help")
        #self.Bind(wx.EVT_MENU, lambda event: self.Destroy(), self.HelpHelp)
        self.Bind(wx.EVT_MENU, self.OnAbout, self.HelpAbout)
        
        #Actually set it
        self.SetMenuBar(self.MenuBar)        
        
    ################################
    #         Event handlers       #
    ################################
    
    def OnOpenConsole(self, event):
        frm = wx.Frame(None)
        from wx.py.crust import Crust
        console = Crust(frm, intro = 'Welcome to the debug console within PDSim', locals = locals())
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(console,1,wx.EXPAND)
        frm.SetSizer(sizer)
        frm.Show()
        
    def OnConfigOpen(self,event):
        
        default_dir = GUIconfig.get('config_path', default = 'configs')
        
        FD = wx.FileDialog(None,"Load Configuration file",defaultDir=default_dir,
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
        FD = wx.FileDialog(None,
                           "Save Configuration file",
                           defaultDir='configs',
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
            
#            for plugin in self.plugins_list:
#                pass
##                if plugin.is_activated():
##                    if hasattr(plugin, ''):
##                        pass
                    
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
        self.Close()
        
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
#            self.MTB.OutputsTB.DataPanel.add_runs([sim])
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
        self.rebuild(module.get_defaults())
        
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