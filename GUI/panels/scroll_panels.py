    # -*- coding: utf-8 -*-

from math import pi, cos, sin
import textwrap

import wx
import wx.grid
from wx.lib.mixins.listctrl import TextEditMixin
from wx.lib.scrolledpanel import ScrolledPanel

import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as WXCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as WXToolbar

import numpy as np

from PDSim.scroll.core import Scroll
from PDSim.scroll.plots import plotScrollSet, ScrollAnimForm
from PDSim.misc.datatypes import AnnotatedValue
import pdsim_panels
from pdsim_panels import LaTeXImageMaker, MotorChoices, PlotPanel
from datatypes import HeaderStaticText, AnnotatedGUIObject
from PDSim.scroll import scroll_geo

# If scipy is available, use its optimization functions, otherwise, 
# use our implementation (for packaging purposes)
try:
    from scipy import optimize
except ImportError:
    import PDSim.misc.solvers as optimize

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
        
class GeometryConverterChoicebook(wx.Choicebook):
    def __init__(self, parent, id=-1, geo = None):
        wx.Choicebook.__init__(self, parent, id)

        self.pagePitch_thickness_height = wx.Panel(self)
        self.AddPage(self.pagePitch_thickness_height,'Pitch, Thickness, Height')

        self.pitch_label = wx.StaticText(self.pagePitch_thickness_height, -1, label = 'Pitch [m]')
        self.thickness_label = wx.StaticText(self.pagePitch_thickness_height, -1, label = 'Thickness [m]')
        self.height_label = wx.StaticText(self.pagePitch_thickness_height, -1, label = 'Height [m]')
        self.W0_label = wx.StaticText(self.pagePitch_thickness_height, -1, label = 'W0 [rad]')
        self.W1_label = wx.StaticText(self.pagePitch_thickness_height, -1, label = 'W1 [rad]')

        self.pitch_value = wx.TextCtrl(self.pagePitch_thickness_height, -1, value = str(geo['pitch']))
        self.thickness_value = wx.TextCtrl(self.pagePitch_thickness_height, -1, value = str(geo['thickness']))
        self.height_value = wx.TextCtrl(self.pagePitch_thickness_height, -1, value = str(geo['height']))
        self.W0_value = wx.TextCtrl(self.pagePitch_thickness_height, -1, value = str(geo['W0']))
        self.W1_value = wx.TextCtrl(self.pagePitch_thickness_height, -1, value = str(geo['W1']))

        sizer_for_outputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)

        # Add all the output objects to the sizer for the outputs
        sizer_for_outputs.AddMany([self.pitch_label, self.pitch_value,
                                   self.thickness_label, self.thickness_value,
                                   self.height_label, self.height_value,
                                   self.W0_label, self.W0_value,
                                   self.W1_label, self.W1_value
                                   ])

        self.pagePitch_thickness_height.SetSizer(sizer_for_outputs)
        sizer_for_outputs.Layout()

    def get_geo(self):
        pitch = float(self.pitch_value.GetValue())
        t = thickness = float(self.thickness_value.GetValue())
        h = height = float(self.height_value.GetValue())
        W1 = float(self.W1_value.GetValue())
        W0 = float(self.W0_value.GetValue())

        rb = base_radius = pitch/(2*pi)
        ro = orbiting_radius = rb*pi - thickness

        # Midline starting wrap angle
        phi_m0 = -W0
        # Initial angles based on offsets off the midline
        phi_i0 = phi_m0+thickness/rb/2.0
        phi_o0 = phi_m0-thickness/rb/2.0
        phi_ie = W1-W0
        phi_is = 0
        phi_os = 0
        
        displacement = -2*pi*h*rb*ro*(3*pi-2*phi_ie+phi_i0+phi_o0)
        volume_ratio = (3*pi-2*phi_ie+phi_i0+phi_o0)/(-2*phi_os-3*pi+phi_i0+phi_o0)

        return dict(displacement = displacement,
                    volume_ratio = volume_ratio,
                    thickness = thickness,
                    orbiting_radius = orbiting_radius,
                    phi_fi0 = phi_i0,
                    phi_fis = phi_is,
                    phi_fos = phi_os,
                    )

class ConvertGeometryFrame(wx.Dialog):
    """ A dialog for converting sets of geometries to the geometry definition used in paper of Bell, IJR, 2013 """
    def __init__(self, geo = None):
        wx.Dialog.__init__(self, None)

        panel = wx.Panel(self)

        self.GCS = GeometryConverterChoicebook(self, geo = geo)
        self.OkButton = wx.Button(self,-1,"Ok")
        self.OkButton.Bind(wx.EVT_BUTTON, self.OnOk)
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.GCS, 1, wx.EXPAND)
        main_sizer.Add(self.OkButton, 1, wx.EXPAND)
        self.SetSizer(main_sizer)
        main_sizer.Layout()

    def OnOk(self, event):
        self.EndModal(wx.ID_OK)

    def get_geo(self):
        """ Get the geometry to be set as a dictionary """
        return self.GCS.get_geo()

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
        
        
        self.CloseButton = wx.Button(panel, label='Close')
        self.CloseButton.Bind(wx.EVT_BUTTON, lambda event: self.Close())
        
        # Do the layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizer_for_outputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.CloseButton, 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        panel.SetSizer(sizer)
        sizer.Layout()
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(panel)
        self.SetSizer(main_sizer)
        main_sizer.Layout()
        self.SetClientSize(main_sizer.GetMinSize())
    
class DischargePortCoordinatesTable(wx.grid.Grid):
    
    def __init__(self, parent, values = None):
        """
        Parameters
        ----------
        parent : wx.window
        values : A 2-element list of lists for all the coordinates (x, y)
        """
        wx.grid.Grid.__init__(self, parent)
        
        #  Make the grid the same shape as the data
        self.CreateGrid(100, 2) # Nrows, Ncolumns
        
        #  Build the headers
        self.SetColLabelValue(0, 'x [m]')
        self.SetColLabelValue(1, 'y [m]')
        
        #  Set the entries in the grid
        if values is not None:
            self.update_from_configfile(values)
    
        #  Bind the events
        self.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.OnCellRightClick)
        
    def OnCellRightClick(self, evt):
        
        # Make a menu
        menu = wx.Menu()
        
        #Build the entries
        menuitem1 = wx.MenuItem(menu, -1, 'Paste from clipboard (Excel format)')
        self.Bind(wx.EVT_MENU, self.OnPaste, menuitem1)
        menu.AppendItem(menuitem1)
        
        if menu.GetMenuItems():
            # Popup the menu.  If an item is selected then its handler
            # will be called before PopupMenu returns.
            self.PopupMenu(menu)
        
        menu.Destroy()
        
    def OnPaste(self, event):
        """
        Paste into the cells in the table
        """
        
        do = wx.TextDataObject()
        if wx.TheClipboard.Open():
            success = wx.TheClipboard.GetData(do)
            wx.TheClipboard.Close()

        data = do.GetText()
        if '\r' in data and '\n' not in data:
            data = data.replace('\r','\n')
        elif '\r\n' in data:
            data = data.replace('\r\n','\n')    
        rows = data.strip().split('\n')
        rows = [row.split('\t') for row in rows]
        
        try:
            for row in rows:
                for el in row:
                    float(el)
            self.update_from_configfile(zip(*rows))
        except ValueError:
            dlg = wx.MessageDialog(None, "Unable to paste from clipboard - bad format")
            dlg.ShowModal()
            dlg.Close()
    
    def ResizeGrid(self, nrows):
        """ Resize the grid to be the right number of rows """
        assert nrows >= 1
        
        if self.GetNumberRows() > nrows:
            while self.GetNumberRows() > nrows:
                self.DeleteRows()
        if self.GetNumberRows() < nrows:
            while self.GetNumberRows() < nrows:
                self.AppendRows()
        
    def update_from_configfile(self, values):
        """
        
        Parameters
        ----------
        values : list of lists, with entries as floating point values
            The first entry is a list (or other iterable) of x values
            
            The second entry is a list (or other iterable) of y values
        """
        self.ResizeGrid(len(values[0]))
                
        for i,(x,y) in enumerate(zip(*values)):
                
            # Values
            self.SetCellValue(i, 0, str(x))
            self.SetCellValue(i, 1, str(y))
                
    def get_coords(self):
        """
        Get the list of lists of values that are used in the table
        """
        x, y = [], []
        for i in range(self.GetNumberRows()):
            x.append(float(self.GetCellValue(i, 0)))
            y.append(float(self.GetCellValue(i, 1)))
        return x, y
    
class DischargePortCoordinatesDialog(wx.Dialog):
    """ A wx.Dialog to hold the grid with the x,y coords """
    def __init__(self, parent, values = None):
        wx.Dialog.__init__(self, parent, title = 'Discharge port coordinates')
        
        self.OKButton = wx.Button(self,label='OK')
        self.OKButton.Bind(wx.EVT_BUTTON, lambda event: self.EndModal(wx.ID_OK))
        self.xy_coords = DischargePortCoordinatesTable(self)
        
        if values is not None:
            self.xy_coords.update_from_configfile(values)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.OKButton, proportion = 0, flag=wx.EXPAND)
        sizer.Add(self.xy_coords, proportion = 1, flag=wx.EXPAND)
        self.SetSizer(sizer)
        sizer.Layout()
        
        w,h = self.GetEffectiveMinSize()
        self.SetSizeWH(w+40,w)
        
        self.xy_coords.ForceRefresh()
        self.Refresh()
        
class DiscCurvesPanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, config):
        pdsim_panels.PDPanel.__init__(self, parent)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        
        if 'disc_curves' in config:
            if 'type' in config['disc_curves']:
                type = config['disc_curves']['type']
            else:
                type = '2Arc'
            if 'r2' in config['disc_curves']:
                r2 = config['disc_curves']['r2']
            else:
                r2 = 0.0
        else:
            type = '2Arc'
            r2 = 0.0

        self.type = wx.Choice(self)
        self.type.AppendItems(['2 Arcs','Arc-Line-Arc'])
        if type == '2Arc':
            self.type.SetSelection(0)
        elif type == 'ArcLineArc':
            self.type.SetSelection(1)
        else:
            raise ValueError
        sizer.Add(self.type)
        
        sizer2.Add(wx.StaticText(self,label='Radius of arc 2'))
        self.r2 = wx.TextCtrl(self, value=str(r2))
        self.r2.SetToolTipString('This can either be the string \"PMP\" for a perfect-meshing profile,\n or alternatively, the second radius r2 in m')
        
        sizer2.Add(self.r2)
        sizer.AddSpacer(3)        
        sizer.Add(sizer2)
        
        self.SetSizer(sizer)
        sizer.Layout()
        
        #  Link callback for refresh of this panel with changing any input
        #  parameter 
        for o in [self.type, self.r2]:
            o.Bind(wx.EVT_KILL_FOCUS, self.GetGrandParent().OnRefresh)
        
geometry_template = """
#  Parameters from the GUI
Vdisp = {Vdisp:s} #[m^3/rev]
Vratio = {Vratio:s} #[-] 
t = {t:s} #[m]
ro = {ro:s} #[m]
phi_i0 = {phi_i0:s} #[rad]
phi_is = {phi_is:s} #[rad]
phi_os = {phi_os:s} #[rad]

#  Set the scroll wrap geometry
sim.set_scroll_geo(Vdisp, # Vdisp [m^3/rev]
               Vratio, # Vratio [-]
               t, # Thickness [m]
               ro, # Orbiting radius [m]
               phi_i0 = phi_i0, # [rad]
               phi_os = phi_os, # [rad]
               phi_is = phi_is) # [rad]
sim.set_disc_geo("{disc_curves_type:s}", r2 = {disc_curves_r2:s})
sim.d_discharge = {d_discharge:s}
{disc_xy_coords_string:s}
sim.geo.delta_flank = {delta_flank:s} # [m]
sim.geo.delta_radial = {delta_radial:s} # [m]

sim.geo.phi_ie_offset = {phi_ie_offset:s}

"""
        
        

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
            
        self.ScrollWrapAnglesButton = wx.Button(scrolled_panel, label = 'View Scroll Wrap Angles')
        self.ScrollWrapAnglesButton.Bind(wx.EVT_BUTTON,self.OnShowWrapGeo)
        self.ConvertGeometryButton = wx.Button(scrolled_panel, label = 'Convert Geometry')
        self.ConvertGeometryButton.Bind(wx.EVT_BUTTON,self.OnConvertGeometry)

        geosizer = wx.BoxSizer(wx.HORIZONTAL)
        geosizer.Add(self.ScrollWrapAnglesButton, 1, wx.ALIGN_CENTER_HORIZONTAL)
        geosizer.Add(self.ConvertGeometryButton, 1, wx.ALIGN_CENTER_HORIZONTAL)

        
        # Build the items and return the list of annotated GUI objects
        annotated_GUI_objects = self.construct_items(annotated_values, 
                                                     sizer = sizer_for_wrap_inputs, 
                                                     parent = scrolled_panel)
                
        #----------------------------------------------------------------------
        # The sizer for all the discharge objects
        sizer_for_discharge_inputs = wx.FlexGridSizer(cols = 1, vgap = 4, hgap = 4)
        
        if 'disc_xy_coords' in config:
            self.disc_xy_coords = config['disc_xy_coords']
            
        self.disc_curves = DiscCurvesPanel(scrolled_panel, config)
        
        self.DiscCoordsButton = wx.Button(scrolled_panel, label = 'Set Disc. Port Coordinates')
        self.DiscCoordsButton.Bind(wx.EVT_BUTTON,self.OnSetDiscPortCoords) 
        self.ClearDiscCoordsButton = wx.Button(scrolled_panel, label = 'Clear Disc. Port Coordinates')
        self.ClearDiscCoordsButton.Bind(wx.EVT_BUTTON,self.OnClearDiscPortCoords) 
        discbut_sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        sizer_for_discharge_inputs.Add(self.disc_curves)
        discbut_sizer.Add(self.DiscCoordsButton, 0, wx.ALIGN_CENTER_HORIZONTAL)
        discbut_sizer.Add(self.ClearDiscCoordsButton, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer_for_discharge_inputs.Add(discbut_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        # Loop over the tube inputs
        annotated_values = []
        for key in ['d_discharge']:
            # Get the annotation and the units for the term 
            annotation, units = self.desc_map[key]
            # Add the annotated object to the list of objects
            annotated_values.append(AnnotatedValue(key, config[key], annotation, units))
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += [self.construct_items(annotated_values,
                                                      sizer = sizer,
                                                      parent = scrolled_panel)]
                                                      
        sizer_for_discharge_inputs.Add(sizer)
        
        #----------------------------------------------------------------------
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
        
        # Link callback for refresh of this panel with changing any input
        # parameter 
        for o in annotated_GUI_objects:
            o.GUI_location.Bind(wx.EVT_KILL_FOCUS, self.OnRefresh)
            
        #Add another callback for use_offset checkbox
        chkbox = self.main.get_GUI_object('use_offset').GUI_location
        chkbox.Bind(wx.EVT_CHECKBOX,self.OnRefresh)
        
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
        sizer.Add(geosizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_wrap_inputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(HeaderStaticText(scrolled_panel, 'Discharge Region Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_discharge_inputs, 0, wx.ALIGN_CENTER_HORIZONTAL)
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
        
    def get_wrap_crossection_involutes(self, axis = 'x'):
        """
        Returns
        phiv : array of centerline involute angles
        
        """
        
        phi0 = (self.Scroll.geo.phi_i0+self.Scroll.geo.phi_o0)/2
        phie = (self.Scroll.geo.phi_ie+self.Scroll.geo.phi_oe)/2
        
        phiv = []
        from PDSim.scroll.scroll_geo import coords_inv
        
        def objective(phi):
            return cos(phi)+(phi-phi0)*sin(phi)
        
        phi = optimize.newton(objective, phi0 + 0.3)
        if phi < phi0: phi += 2*pi
        
        while phi < phie:
            phiv.append(phi)
            phi = optimize.newton(objective, phi + pi)
        
        phiv.append(phi)
        
        return phiv, self.Scroll.geo.h, self.Scroll.geo.t
        
    def OnSetDiscPortCoords(self, event = None):
        
        #  Get the current values for the discharge port coordinates
        if hasattr(self,'disc_xy_coords'):
            values = self.disc_xy_coords
        else:
            values = None
            
        dlg = DischargePortCoordinatesDialog(None, values)
        if dlg.ShowModal() == wx.ID_OK:
            x,y = dlg.xy_coords.get_coords()
            self.disc_xy_coords = x,y
            self.OnRefresh()
        dlg.Destroy()
    
    def OnClearDiscPortCoords(self, event = None):
        """ Clear the coordinates for the discharge port """
        del self.disc_xy_coords
        self.OnRefresh()
        
    def OnShowWrapGeo(self, event = None):
        if event is not None: event.Skip()
        
        frm = ScrollWrapAnglesFrame(self.Scroll.geo)
        frm.Show()

    def OnConvertGeometry(self, event = None):
        if event is not None: event.Skip()
    
        def get(key):
            # Compact code to get a parameter from the main database
            return self.main.get_GUI_object_value(key)

        W0 = self.Scroll.geo.t/(2*self.Scroll.geo.rb)-self.Scroll.geo.phi_fi0
        W1_minus_W0 = self.Scroll.geo.phi_fie-self.Scroll.geo.phi_fis
        geo = dict(pitch = 2*pi*self.Scroll.geo.rb,
                   thickness = self.Scroll.geo.t,
                   height = self.Scroll.geo.h,
                   W0 = W0,
                   W1 = W1_minus_W0 + W0
                   )

        frm = ConvertGeometryFrame(geo = geo)
        if frm.ShowModal() == wx.ID_OK:
            geo = frm.get_geo()

            self.main.get_GUI_object('t').SetValue(str(geo['thickness']))
            self.main.get_GUI_object('Vdisp').SetValue(str(geo['displacement']))
            self.main.get_GUI_object('Vratio').SetValue(str(geo['volume_ratio']))
            self.main.get_GUI_object('ro').SetValue(str(geo['orbiting_radius']))
            for key in ['phi_fi0','phi_fos','phi_fis']:
                self.main.get_GUI_object(key).SetValue(str(geo[key]))
        frm.Destroy()
        
    def OnAnimate(self, event = None):
        
        if hasattr(self,'disc_xy_coords'):
            disc_xy_coords = self.disc_xy_coords
        else:
            disc_xy_coords = None 
            
        SAF = ScrollAnimForm(self.Scroll.geo, size=(400,400), param_dict = self.main.get_GUI_object_value_dict(), disc_xy_coords = disc_xy_coords)
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
        
        if self.disc_curves.type.GetStringSelection() == '2 Arcs':
            disc_curves_type = '2Arc'
        elif self.disc_curves.type.GetStringSelection() == 'Arc-Line-Arc':
            disc_curves_type = 'ArcLineArc'
        else:
            raise ValueError
        
        #  Get r2 as a string, convert to a floating point value if possible
        r2 = self.disc_curves.r2.GetValue()
        try:
            r2 = float(r2)
        except ValueError:
            pass 
        
        self.Scroll.set_disc_geo(disc_curves_type, r2 = r2)
        
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
            if not hasattr(self,'disc_xy_coords'):
                t = np.linspace(0, 2*np.pi)
                x = self.Scroll.geo.xa_arc1 + d_discharge/2*np.cos(t)
                y = self.Scroll.geo.ya_arc1 + d_discharge/2*np.sin(t)
                self.ax.plot(x,y,'--')
        except KeyError:
            pass
        
        if hasattr(self,'disc_xy_coords'):
            self.ax.plot(self.disc_xy_coords[0],self.disc_xy_coords[1]) 
            
        self.PP.canvas.draw()
        
    def get_config_chunk(self):
        
        #  All the conventional terms
        keys = ['Vdisp','Vratio','t','ro','phi_fi0','phi_fis','phi_fos',
                'use_offset','delta_offset','delta_flank','delta_radial',
                'd_discharge','inlet_tube_length', 'inlet_tube_ID', 
                'outlet_tube_length', 'outlet_tube_ID']
        
        #  Dictionary of the values
        d = {key:self.main.get_GUI_object_value(key) for key in keys}
        
        #  Added values for the discharge curves
        disc_type = self.disc_curves.type.GetStringSelection()
        if disc_type == '2 Arcs':
            disc_type = '2Arc'
        elif disc_type == 'Arc-Line-Arc':
            disc_type = 'ArcLineArc'
        else:
            raise ValueError
        disc_r2 = self.disc_curves.r2.GetValue()
        d.update(dict(disc_curves = dict(type = disc_type, r2 = disc_r2)))
        
        #  Added values for the discharge port curves
        if hasattr(self,'disc_xy_coords'):
            d.update(dict(disc_xy_coords = self.disc_xy_coords))
        return d
        
    def get_script_chunks(self, plugin_chunks = None):
        
        def get(key):
            # Compact code to get a parameter from the main database
            return self.main.get_GUI_object(key).GetValue()
    
        if get('use_offset'):
            phi_ie_offset = str(pi)
        else:
            phi_ie_offset = str(0)
            
        template = "sim.geo.xvec_disc_port = np.array({x:s})\nsim.geo.yvec_disc_port = np.array({y:s})"
        if hasattr(self, 'disc_xy_coords'): 
            disc_xy_coords_string = textwrap.dedent(template.format(x = str(self.disc_xy_coords[0]),
                                                                    y = str(self.disc_xy_coords[1]))
                                                    )
        else:
            disc_xy_coords_string = ''
            
        if self.disc_curves.type.GetStringSelection() == '2 Arcs':
            disc_curves_type = '2Arc'
        elif self.disc_curves.type.GetStringSelection() == 'Arc-Line-Arc':
            disc_curves_type = 'ArcLineArc'
        else:
            raise ValueError
        
        #  Get r2 as a string, convert to a floating point value if possible
        r2 = self.disc_curves.r2.GetValue()
        try:
            #  If this works, r2 is a floating point expressed as a string,
            #  leave it alone
            float(r2)
        except ValueError:
            #  r2 is PMP, wrap it in quotes
            r2 = '"' + r2 + '"'
        
        #Parameters to be set in the string:
        str_params = dict(Vdisp = get('Vdisp'),
                          Vratio = get('Vratio'),
                          t = get('t'),
                          ro = get('ro'),
                          phi_i0 = get('phi_fi0'),
                          phi_os = get('phi_fos'),
                          phi_is = get('phi_fis'),
                          delta_flank = get('delta_flank'),
                          delta_radial = get('delta_radial'),
                          d_discharge = get('d_discharge'),
                          disc_curves_type = disc_curves_type,
                          disc_curves_r2 = r2,
                          phi_ie_offset = phi_ie_offset,
                          disc_xy_coords_string = disc_xy_coords_string
                          )

        core = textwrap.dedent(geometry_template.format(**str_params))
        
        # Add plugin injected chunks
        allowed = ['ScrollGeometryPanel_After', 'ScrollGeometryPanel_Before']
        if isinstance(plugin_chunks,dict):
            for key,chunk in plugin_chunks.iteritems():
                if key in allowed:
                    core += chunk
                    
        return core
        
        
class FlowOptions(pdsim_panels.PDPanel):
    """
    Takes a list of dictionaries in and creates a panel with a dropdown to select
    the model and a set of objects to change the parameters
    
    Returns
    -------
    A list of annotated GUI objects for each item that is created
    """
    def __init__(self, parent, pathname, choices_list, register_objects = True):
        wx.Panel.__init__(self, parent)
        
        annotated_objects = []
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        self.choice_book = wx.Choicebook(self, -1)
        
        for choice in choices_list:
            panel = wx.Panel(self.choice_book)
            self.choice_book.AddPage(panel, text = choice['model'])
            
            panel_sizer = wx.FlexGridSizer(cols = 2)
            panel_annotated_objects = []
            
            for option in choice['options']:
                term_name = 'flow path' + pathname + '|' + option
                value = choice['options'][option]
                panel_annotated_objects.append(AnnotatedValue(term_name, value, term_name, ''))
            
            #Annotated GUI objects
            panel_AGO = self.construct_items(panel_annotated_objects,
                                             sizer = panel_sizer,
                                             parent = panel)
            
            if register_objects:
                self.GetTopLevelParent().register_GUI_objects(panel_AGO)
                
            panel.SetSizer(panel_sizer)
            panel_sizer.Layout()
            
        sizer.Add(self.choice_book, 0)
        self.SetSizer(sizer)
        sizer.Layout()
                
class MassFlowPanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, configdict, **kwargs):
    
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        options = {}
        for flow in ['sa-s1', 'sa-s2', 'inlet.2-sa','d1-dd','d2-dd']:
            if flow not in configdict:
                options[flow] = dict(model = 'IsentropicNozzle', options = dict(Xd = 0.7))
            else:
                options[flow] = dict(model = configdict[flow]['model'], options = configdict[flow]['options'])
            
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.flow1 = FlowOptions(self, 'sa-s1', [options['sa-s1']])
        self.flow2 = FlowOptions(self, 'sa-s2', [options['sa-s2']])
        self.flow3 = FlowOptions(self, 'inlet.2-sa', [options['inlet.2-sa']])
        self.flow4 = FlowOptions(self, 'd1-dd', [options['d1-dd']])
        self.flow5 = FlowOptions(self, 'd2-dd', [options['d2-dd']])
        
        sizer.Add(pdsim_panels.HeaderStaticText(self,'Flow model parameters') , 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(10)
        sizer.Add(wx.StaticText(self,label='sa-s1'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.flow1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(wx.StaticText(self,label='sa-s2'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.flow2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(wx.StaticText(self,label='inlet.2-sa'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.flow3, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(wx.StaticText(self,label='d1-dd'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.flow4, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(wx.StaticText(self,label='d2-dd'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.flow5, 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        self.SetSizer(sizer)
        sizer.Layout()
    
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
        
    def get_config_chunk(self):
        
        Xd_sa_s1 = self.main.get_GUI_object_value('flow pathsa-s1|Xd')
        Xd_sa_s2 = self.main.get_GUI_object_value('flow pathsa-s2|Xd')
        Xd_inlet = self.main.get_GUI_object_value('flow pathinlet.2-sa|Xd')
        Xd_d1_dd = self.main.get_GUI_object_value('flow pathd1-dd|Xd')
        Xd_d2_dd = self.main.get_GUI_object_value('flow pathd2-dd|Xd')
                       
        configdict = {}
        configdict['sa-s1'] = dict(options = dict(Xd = Xd_sa_s1), 
                                   model='IsentropicNozzle')
        configdict['sa-s2'] = dict(options = dict(Xd = Xd_sa_s2), 
                                   model='IsentropicNozzle')
        configdict['inlet.2-sa'] = dict(options = dict(Xd = Xd_inlet), 
                                        model='IsentropicNozzle')
        configdict['d1-dd'] = dict(options = dict(Xd = Xd_d1_dd), 
                                        model='IsentropicNozzle')
        configdict['d2-dd'] = dict(options = dict(Xd = Xd_d2_dd), 
                                        model='IsentropicNozzle')
        
        return configdict
        
    def get_script_chunks(self, plugin_chunks = None):
        
        Xd_dict = dict(Xd_sa_s1 = str(self.main.get_GUI_object_value('flow pathsa-s1|Xd')),
                       Xd_sa_s2 = str(self.main.get_GUI_object_value('flow pathsa-s2|Xd')),
                       Xd_inlet = str(self.main.get_GUI_object_value('flow pathinlet.2-sa|Xd')),
                       Xd_d1_dd = str(self.main.get_GUI_object_value('flow pathd1-dd|Xd')),
                       Xd_d2_dd = str(self.main.get_GUI_object_value('flow pathd2-dd|Xd')),
                       inlet_tube_length = str(self.main.get_GUI_object_value('inlet_tube_length')),
                       outlet_tube_length = str(self.main.get_GUI_object_value('outlet_tube_length')),
                       inlet_tube_ID = str(self.main.get_GUI_object_value('inlet_tube_ID')),
                       outlet_tube_ID = str(self.main.get_GUI_object_value('outlet_tube_ID')),
                       )
          
        return textwrap.dedent(
            """
            # Add all the control volumes
            sim.auto_add_CVs(inletState, outletState)
            
            # Get the guess for the mass flow rate
            mdot_guess = inletState.rho*sim.Vdisp*sim.omega/(2*pi)
    
            # Add both the inlet and outlet tubes
            sim.add_tube(Tube(key1 = 'inlet.1',
                              key2 = 'inlet.2',
                              L = {inlet_tube_length:s},
                              ID = {inlet_tube_ID:s},
                              mdot = mdot_guess, 
                              State1 = inletState.copy(),
                              fixed = 1,
                              TubeFcn = sim.TubeCode))
            sim.add_tube(Tube(key1 = 'outlet.1',
                              key2 = 'outlet.2',
                              L = {outlet_tube_length:s},
                              ID = {outlet_tube_ID:s},
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
            FP.A = pi*{inlet_tube_ID:s}**2/4*{Xd_inlet:s}
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
                        
            sim.add_flow(FlowPath(key1 = 'outlet.1',
                                 key2 = 'dd',
                                 MdotFcn = sim.DISC_DD,
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                                 )
                        )
        
            sim.add_flow(FlowPath(key1 = 'outlet.1',
                                 key2 = 'ddd',
                                 MdotFcn = sim.DISC_DD,
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                               )
                        )
            sim.add_flow(FlowPath(key1 = 'outlet.1',
                                 key2 = 'd1',
                                 MdotFcn = sim.DISC_D1,
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                                 )
                        )
            sim.add_flow(FlowPath(key1 = 'outlet.1',
                                 key2 = 'ddd',
                                 MdotFcn = sim.DISC_D1,
                                 MdotFcn_kwargs = dict(X_d = 0.7)
                                 )
                        )
            
            sim.add_flow(FlowPath(key1='d1',
                                  key2='dd',
                                  MdotFcn=sim.D_to_DD,
                                  MdotFcn_kwargs = dict(X_d = {Xd_d1_dd:s})
                                  )
                              )
            sim.add_flow(FlowPath(key1='d2',
                                  key2='dd',
                                  MdotFcn=sim.D_to_DD,
                                  MdotFcn_kwargs = dict(X_d = {Xd_d2_dd:s})
                                  )
                        )
            """.format(**Xd_dict)
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
    
class OSCrossSectionFrame(wx.Frame):
    def __init__(self, dictionary, phiv, h, w):
        """
        Parameters
        ----------
        dictionary : dict
            Dictionary from the GUI of all the annotated terms
        """
        wx.Frame.__init__(self,None)
        
        from PDSim.scroll.plots import OSCrossSectionPanel
        panel = OSCrossSectionPanel(self, dictionary, phiv, h, w)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(panel, 0, wx.EXPAND)
        self.SetSizer(sizer)
        sizer.Layout()
        
        self.SetSize(sizer.GetMinSize())

class MechanicalLossesChoices(wx.Choicebook):
    def __init__(self, parent):
        wx.Choicebook.__init__(self, parent, -1)
        
        self.page_mech_normal=wx.Panel(self)
        self.AddPage(self.page_mech_normal,'Normal Mechanical Losses')
        
        self.page_spec_eta_mech=wx.Panel(self)
        self.AddPage(self.page_spec_eta_mech,'Specified Mechanical Efficiency')
        
        self.page_spec_mech_losses=wx.Panel(self)
        self.AddPage(self.page_spec_mech_losses, 'Specified Mechanical Losses')
        
        
class MechanicalLossesPanel(pdsim_panels.PDPanel):
    
    desc_map = dict(h_shell = ('Shell-ambient mean HTC [kW/m\xb2/K]','kW/m^2/K'),
                    A_shell = ('Shell outer area [m\xb2]','m^2'),
                    Tamb = ('Ambient temperature [K]','K'),
                    mu_oil = ('Viscosity of the oil [Pa-s]','Pa-s'),
                    D_upper_bearing = ('Upper bearing journal diameter [m]','m'),
                    L_upper_bearing = ('Upper bearing length [m]','m'),
                    c_upper_bearing = ('Upper bearing clearance [m]','m'),
                    D_crank_bearing = ('Crank bearing journal diameter [m]','m'),
                    L_crank_bearing = ('Crank bearing length [m]','m'),
                    c_crank_bearing = ('Crank bearing clearance [m]','m'),
                    D_lower_bearing = ('Lower bearing journal diameter [m]','m'),
                    L_lower_bearing = ('Lower bearing length [m]','m'),
                    c_lower_bearing = ('Lower bearing clearance [m]','m'),
                    journal_tune_factor = ('Tuning factor on journal bearing losses [-]','-'),
                    thrust_friction_coefficient = ('Thrust bearing friction coefficient [-]','-'),
                    thrust_ID = ('Thrust bearing inner diameter [m]','m'),
                    thrust_OD = ('Thrust bearing outer diameter [m]','m'),  
                    L_ratio_bearings = ('Ratio of lengths to the bearings [-]','-'),
                    scroll_plate_thickness = ('Thickness of the orbiting scroll plate [m]','m',0.002),
                    scroll_plate_diameter = ('Effective diameter of the orbiting scroll plate [m]','m',0.014),
                    scroll_density = ('Orbiting scroll material density [kg/m\xb3]','kg/m^3',2700),
                    scroll_added_mass = ('Additional OS mass added at COM [kg]','kg',0.0),
                    oldham_ring_radius = ('Oldham ring radius [m]','m',0.06),
                    oldham_mass = ('Mass of the Oldham ring [kg]','kg',0.1),
                    oldham_thickness = ('Height of the Oldham ring (without the keys) [m]','m',0.008),
                    oldham_key_height = ('Height of the keys of the Oldham ring [m]','m',0.006),
                    oldham_key_width = ('Width of the keys of the Oldham ring [m]','m',0.006),
                    oldham_key_friction_coefficient = ('Friction coefficient of the Oldham ring [-]','-',0.01),
                    oldham_rotation_beta = ('Angle between Oldham sliding axis and x-axis [radian]','rad',0),
                    HTC = ('Heat transfer coefficient in the scrolls [kW/m\xb2/K]','kW/m^2/K'),
                    detailed_analysis = ('Use detailed analysis of the mechanical losses','',True),
                    suction_fraction = ('Fraction of motor losses to suction gas','',1.0),
                    pin1_ybeta_offset = ('Offset of pin #1 in +y_beta direction [m]','m',0.0),
                    pin2_ybeta_offset = ('Offset of pin #2 in +y_beta direction [m]','m',0.0),
                    pin3_xbeta_offset = ('Offset of pin #3 in +x_beta direction [m]','m',0.0),
                    pin4_xbeta_offset = ('Offset of pin #4 in +x_beta direction [m]','m',0.0),
                    specified_mechanical_efficiency = ('Specified mechanical efficiency [-]','-',0.9),
                    specified_mechanical_losses_kW = ('Specified mechanical losses [kW]','kW',0.0)
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
        
        if 'orbiting_scroll_mass' in config:
            import warnings
            warnings.warn('the term "orbiting_scroll_mass" has been deprecated, please remove it from your configuration')
            config.pop('orbiting_scroll_mass')
        
        self.motor_choices = MotorChoices(scrolled_panel)
        
        if ('eta_motor' in config
            and 'eta_motor_coeffs' not in config
            and 'tau_motor_coeffs' not in config
            and 'omega_motor_coeffs' not in config):
            eta_motor = config['eta_motor']
            #Only eta_motor is provided, use it in the motor panel
            self.motor_choices.SetSelection(0)
            #Set the value in the panel
            self.motor_choices.eta_motor.SetValue(str(eta_motor))
            # When the motor efficiency is changed by something else, it means
            # we want to use the motor efficiency rather than the motor curves,
            # so set it back to using constant efficiency
            self.motor_choices.eta_motor.Bind(wx.EVT_TEXT,lambda event: self.motor_choices.SetSelection(0))
            
            AGO_motor = AnnotatedGUIObject(AnnotatedValue('eta_motor', eta_motor, 'Motor Efficiency [-]','-'),self.motor_choices.eta_motor)
            
            self.main.register_GUI_objects(AGO_motor)
            
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
            
        self.keys_for_config = []
        
        self.mechanical_model_choices = MechanicalLossesChoices(scrolled_panel)
        
        if 'specified_mechanical_efficiency' in config:
            self.mechanical_model_choices.SetSelection(1)
        elif 'specified_mechanical_losses_kW' in config:
            self.mechanical_model_choices.SetSelection(2)
        
        #----------------------------------------------------------------------
        # The sizer for all the specified mechanical efficiency terms
        sizer_for_spec_etamech_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the inputs
        annotated_values = self.get_annotated_values(['specified_mechanical_efficiency'])
        
        # Build the items and return the list of annotated GUI objects, add to existing list
        AGO = self.construct_items(annotated_values,
                                   sizer = sizer_for_spec_etamech_inputs,
                                   parent = self.mechanical_model_choices.page_spec_eta_mech)
        self.main.register_GUI_objects(AGO)
                                                      
        #----------------------------------------------------------------------
        # The sizer for all the specified mechanical losses terms
        sizer_for_spec_mech_losses_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the inputs
        annotated_values = self.get_annotated_values(['specified_mechanical_losses_kW'])
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        AGO = self.construct_items(annotated_values,
                                   sizer = sizer_for_spec_mech_losses_inputs,
                                   parent = self.mechanical_model_choices.page_spec_mech_losses)
        self.main.register_GUI_objects(AGO)
                                                      
        #----------------------------------------------------------------------
        # The sizer for all the orbiting scroll terms
        sizer_for_orbiting_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the inputs
        keys = ['scroll_plate_thickness', 'scroll_plate_diameter','scroll_density', 'scroll_added_mass']
        annotated_values = self.get_annotated_values(keys)
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects = self.construct_items(annotated_values,
                                                     sizer = sizer_for_orbiting_inputs,
                                                     parent = scrolled_panel)
        self.main.register_GUI_objects(annotated_GUI_objects)
        
        self.MassButton = wx.Button(scrolled_panel,label='Calculate')
        sizer_for_orbiting_inputs.Add(wx.StaticText(scrolled_panel,label = 'Orbiting Scroll Mass [kg]'))
        sizer_for_orbiting_inputs.Add(self.MassButton)
        self.MassButton.Bind(wx.EVT_BUTTON,self.OnCalculateScrollMass)
        
        annotated_GUI_objects = []
        self.config = config
        
        
        #----------------------------------------------------------------------
        # The sizer for all the heat transfer terms
        sizer_for_HT_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the HT inputs
        annotated_values = self.get_annotated_values(['h_shell','A_shell','Tamb','HTC','suction_fraction'])
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                       sizer = sizer_for_HT_inputs,
                                                       parent = scrolled_panel)
        
        #----------------------------------------------------------------------
        # The sizer for all the journal bearings terms
        sizer_for_journal_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        keys = ['D_upper_bearing','L_upper_bearing','c_upper_bearing',
                'D_crank_bearing','L_crank_bearing','c_crank_bearing',
                'D_lower_bearing','L_lower_bearing','c_lower_bearing',
                'journal_tune_factor','L_ratio_bearings']
        annotated_values = self.get_annotated_values(keys)
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                      sizer = sizer_for_journal_inputs,
                                                      parent = self.mechanical_model_choices.page_mech_normal)
        
        #----------------------------------------------------------------------
        # The sizer for all the Oldham ring terms
        sizer_for_oldham_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the oldham inputs

        keys = ['oldham_mass', 'oldham_thickness', 'oldham_key_height', 
                'oldham_key_width', 'oldham_key_friction_coefficient',
                'oldham_rotation_beta','oldham_ring_radius','pin1_ybeta_offset',
                'pin2_ybeta_offset','pin3_xbeta_offset','pin4_xbeta_offset']
        annotated_values = self.get_annotated_values(keys)
        
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                      sizer = sizer_for_oldham_inputs,
                                                      parent = self.mechanical_model_choices.page_mech_normal)
        
        #----------------------------------------------------------------------
        # The sizer for all the thrust bearing terms
        sizer_for_thrust_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the inputs
        keys = ['thrust_friction_coefficient', 'thrust_ID', 'thrust_OD']
        annotated_values = self.get_annotated_values(keys)
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                      sizer = sizer_for_thrust_inputs,
                                                      parent = self.mechanical_model_choices.page_mech_normal)
        
        #----------------------------------------------------------------------
        # The sizer for all the general bearing terms
        sizer_for_general_inputs = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        # Loop over the inputs
        annotated_values = self.get_annotated_values(['mu_oil','detailed_analysis'])
            
        # Build the items and return the list of annotated GUI objects, add to existing list
        annotated_GUI_objects += self.construct_items(annotated_values,
                                                      sizer = sizer_for_general_inputs,
                                                      parent = self.mechanical_model_choices.page_mech_normal)
                                                      
        
        
        # Register terms in the GUI database
        self.main.register_GUI_objects(annotated_GUI_objects)
        
        self.main.get_GUI_object('L_ratio_bearings').GUI_location.SetToolTipString('Ratio of z1/z2, where\n\nz1 : the length from the centerline of the upper bearing to the lower bearing\nz2 : the length from the centerline of the upper bearing to the orbiting scroll bearing')
        
        self.ViewButton = wx.Button(scrolled_panel, label='View Cross-Section')
        self.ViewButton.Bind(wx.EVT_BUTTON, self.OnViewCrossSection)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(HeaderStaticText(self.mechanical_model_choices.page_mech_normal, 'General Mechanical Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_general_inputs,0,wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        sizer.Add(HeaderStaticText(self.mechanical_model_choices.page_mech_normal, 'Bearing Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_journal_inputs,0,wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        sizer.Add(HeaderStaticText(self.mechanical_model_choices.page_mech_normal, 'Oldham Ring Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_oldham_inputs,0,wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        sizer.Add(HeaderStaticText(self.mechanical_model_choices.page_mech_normal, 'Thrust Bearing Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_thrust_inputs,0,wx.ALIGN_CENTER_HORIZONTAL)
        self.mechanical_model_choices.page_mech_normal.SetSizer(sizer)
        sizer.Layout()
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizer_for_spec_etamech_inputs,0,wx.ALIGN_CENTER_HORIZONTAL)
        self.mechanical_model_choices.page_spec_eta_mech.SetSizer(sizer)
        sizer.Layout()
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizer_for_spec_mech_losses_inputs,0,wx.ALIGN_CENTER_HORIZONTAL)
        self.mechanical_model_choices.page_spec_mech_losses.SetSizer(sizer)
        sizer.Layout()
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.ViewButton,0,wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(10)
        sizer.Add(HeaderStaticText(scrolled_panel, "Motor Model"), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.motor_choices,0,wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        sizer.Add(HeaderStaticText(scrolled_panel, 'Heat Transfer Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_HT_inputs,0,wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        sizer.Add(HeaderStaticText(scrolled_panel, 'Orbiting Scroll Inputs'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(sizer_for_orbiting_inputs,0,wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        sizer.Add(HeaderStaticText(scrolled_panel, 'Mechanical Loss Models'), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(5)
        sizer.Add(self.mechanical_model_choices,0,wx.ALIGN_CENTER_HORIZONTAL)
        sizer.AddSpacer(20)
        scrolled_panel.SetSizer(sizer)
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)
        self.SetSizer(main_sizer)
        sizer.Layout()
                    
    def get_config_chunk(self):
        configdict = {}
        keys_for_config = list(self.keys_for_config)
        
        if self.motor_choices.GetSelection() == 0:
            configdict['eta_motor'] = float(self.motor_choices.eta_motor.GetValue())
        elif self.motor_choices.GetSelection() == 1:
            c = self.motor_choices.MCT.get_coeffs()
            configdict['tau_motor_coeffs'] = c[0]
            configdict['eta_motor_coeffs'] = c[1]
            configdict['omega_motor_coeffs'] = c[2]
            
        if self.mechanical_model_choices.GetSelection() == 0:
            for key in ['specified_mechanical_efficiency',
                        'specified_mechanical_losses_kW']:
                if key in keys_for_config:
                    keys_for_config.pop(keys_for_config.index(key))
        elif self.mechanical_model_choices.GetSelection() == 1:
            for key in ['specified_mechanical_losses_kW']:
                if key in keys_for_config:
                    keys_for_config.pop(keys_for_config.index(key))
        elif self.mechanical_model_choices.GetSelection() == 2:
            for key in ['specified_mechanical_efficiency']:
                if key in keys_for_config:
                    keys_for_config.pop(keys_for_config.index(key))
            
        for key in keys_for_config:
            configdict[key] = self.main.get_GUI_object_value(key)
            
        return configdict
                           
    def get_script_chunks(self, plugin_chunks = None):
        """
        Returns a formatted string for the script that will be execfile-d
        """

        if self.motor_choices.GetSelection() == 0:
            #Use the value for the motor efficiency
            motor_chunk = textwrap.dedent(
               """
               sim.motor = Motor()
               sim.motor.set_eta({eta_motor:s})
               sim.motor.suction_fraction = {suction_fraction:g}
               
               from PDSim.core.core import struct
               sim.mech = struct()
               """.format(eta_motor = self.motor_choices.eta_motor.GetValue(),
                          suction_fraction = self.main.get_GUI_object_value('suction_fraction')
                          )
                                          )
        elif self.motor_choices.GetSelection() == 1:
            # Get the tuple of list of coeffs from the MCT, then unpack the tuple
            # back into the call to set the coefficients
            c = self.motor_choices.MCT.get_coeffs()
            #Will set the type flag itself
            motor_chunk = textwrap.dedent(
               """
               sim.motor = Motor()
               sim.motor.set_coeffs(tau_coeffs = {tau_coeffs:s},
                                    eta_coeffs = {eta_coeffs:s},
                                    omega_coeffs = {omega_coeffs:s})
               sim.motor.suction_fraction = {suction_fraction:g}
               
               from PDSim.core.core import struct
               sim.mech = struct()
               """.format(tau_coeffs = str(c[0]),
                          eta_coeffs = str(c[1]),
                          omega_coeffs = str(c[2]),
                          suction_fraction = self.main.get_GUI_object_value('suction_fraction')
                          )        
                                          )
        else:
            raise NotImplementedError
        
        #Terms that do not go in the mech struct
        for term in ['h_shell','A_shell','Tamb','HTC']:
            val = self.main.get_GUI_object_value(term)
            motor_chunk += 'sim.{name:s} = {value:s}\n'.format(name = term,
                                                               value = str(val))
        
        # Terms that always go in the mech struct
        for term in ['scroll_plate_thickness','scroll_plate_diameter','scroll_added_mass','scroll_density']:
            val = self.main.get_GUI_object_value(term)
            motor_chunk += 'sim.mech.{name:s} = {value:s}\n'.format(name = term,
                                                                    value = str(val))
        
        if self.mechanical_model_choices.GetSelection() == 0:
            #Terms that go in the mech struct
            for term in ['mu_oil','detailed_analysis','journal_tune_factor',
                    'D_upper_bearing','L_upper_bearing','c_upper_bearing',
                    'D_crank_bearing','L_crank_bearing','c_crank_bearing',
                    'D_lower_bearing','L_lower_bearing','c_lower_bearing',
                    'thrust_friction_coefficient', 'thrust_ID', 'thrust_OD', 
                    'L_ratio_bearings', 'oldham_key_friction_coefficient', 
                    'oldham_ring_radius', 'oldham_key_width', 'oldham_mass', 
                    'oldham_thickness', 'oldham_key_height','oldham_rotation_beta',
                    'pin1_ybeta_offset','pin2_ybeta_offset','pin3_xbeta_offset','pin4_xbeta_offset'
                    ]:
                val = self.main.get_GUI_object_value(term)
                motor_chunk += 'sim.mech.{name:s} = {value:s}\n'.format(name = term,
                                                                        value = str(val))
        elif self.mechanical_model_choices.GetSelection() == 1:
            for term in ['specified_mechanical_efficiency']:
                val = self.main.get_GUI_object_value(term)
                motor_chunk += 'sim.mech.{name:s} = {value:s}\n'.format(name = term,
                                                                        value = str(val)) 
        elif self.mechanical_model_choices.GetSelection() == 2:
            for term in ['specified_mechanical_losses_kW']:
                val = self.main.get_GUI_object_value(term)
                motor_chunk += 'sim.mech.{name:s} = {value:s}\n'.format(name = term,
                                                                        value = str(val)) 
        
        # Handle the orbiting scroll mass plus any additional mass
        motor_chunk += 'm, zcm = sim.calculate_scroll_mass()\nsim.mech.orbiting_scroll_mass = m\nsim.mech.scroll_zcm__thrust_surface = zcm\n'
        
        return motor_chunk
    
    def OnViewCrossSection(self, event):
        
        # Get the panel that has the geometry parameters 
        GeoPanel = self.main.get_GUI_object('Vratio').GUI_location.GetGrandParent()
        
        phiv, h, w = GeoPanel.get_wrap_crossection_involutes()
        frm = OSCrossSectionFrame(self.main.get_GUI_object_value_dict(), phiv, h, w)
        frm.Show()
        
    def OnCalculateScrollMass(self, event):
        #mtotal,zcm = self.calculate_scroll_mass()
        #template = """Scroll Mass : {m:g} kg\nCentroid : {zcm:g} m (relative to the thrust surface)"""
        #dlg = wx.MessageDialog(None, template.format(m=mtotal, zcm = zcm))
        dlg = wx.MessageDialog(None, 'Temporarily disabled - sorry')
        dlg.ShowModal()
        dlg.Destroy()
        
class InvoluteToCoords(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, title = 'Involute to coordinates')
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        FGS = wx.FlexGridSizer(cols = 2, vgap = 4, hgap = 4)
        
        self.inv = wx.Choice(self)
        self.inv.AppendItems(['Fixed Inner','Fixed Outer'])
        self.angle = wx.TextCtrl(self)
        self.offset = wx.TextCtrl(self)
        self.AddButton = wx.Button(self, label = 'Add')
        self.AddButton.Bind(wx.EVT_BUTTON, self.OnAdd)
        
        FGS.Add(wx.StaticText(self,label='Involute'))
        FGS.Add(self.inv)
        FGS.Add(wx.StaticText(self,label='Angle [rad]'))
        FGS.Add(self.angle)
        FGS.Add(wx.StaticText(self,label='Offset [m]'))
        FGS.Add(self.offset)
        
        sizer.Add(FGS)
        sizer.Add(self.AddButton)
        
        self.SetSizer(sizer)
        self.Fit()
        
        #Bind a key-press event to all objects to get Esc 
        children = self.GetChildren()
        for child in children:
            child.Bind(wx.EVT_KEY_UP,  self.OnKeyPress)
    
    def OnAdd(self, event = None):
        self.EndModal(wx.ID_OK)
        
    def OnKeyPress(self,event = None):
        """ cancel if Escape key is pressed """
        event.Skip()
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self.EndModal(wx.ID_CANCEL)
        
class AddSensorDialog(wx.Dialog):
    
    def __init__(self, parent, geo):
        
        wx.Dialog.__init__(self, parent, title = 'Virtual Sensor Selection', size = (500,500))
        
        # local copy of geometry
        self.geo = geo
        
        # The plot of the scroll wraps
        self.PP = PlotPanel(self)
        self.ax = self.PP.figure.add_axes((0, 0, 1, 1))
        
        self.FromInvolute = wx.Button(self, label = 'From Involute...')
        self.FromInvolute.Bind(wx.EVT_BUTTON, self.OnFromInvolute)
        self.Accept = wx.Button(self, label = 'Accept')
        self.Accept.Bind(wx.EVT_BUTTON, self.OnAccept)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        xsizer = wx.BoxSizer(wx.HORIZONTAL)
        ysizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.x = wx.TextCtrl(self, value='')
        self.y = wx.TextCtrl(self, value='')
        
        xsizer.Add(wx.StaticText(self,label='x [m]'))
        xsizer.Add(self.x)
        ysizer.Add(wx.StaticText(self,label='y [m]'))
        ysizer.Add(self.y)
        
        sizer.Add(self.PP, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(xsizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(ysizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.FromInvolute, 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.Accept, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.SetSizer(sizer)
        
        self.OnRefresh()
        self.Fit()
        
        #Bind a key-press event to all objects to get Esc 
        children = self.GetChildren()
        for child in children:
            child.Bind(wx.EVT_KEY_UP,  self.OnKeyPress)
    
    def OnAccept(self, event = None):
        self.EndModal(wx.ID_OK)
        
    def OnKeyPress(self,event = None):
        """ cancel if Escape key is pressed """
        event.Skip()
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self.EndModal(wx.ID_CANCEL)
        
    def OnFromInvolute(self, event = None):
        key_dict = {'Orbiting Inner': 'oi', 'Orbiting Outer':'oo','Fixed Inner':'fi','Fixed Outer':'fo'}
        dlg = InvoluteToCoords(None)
        if dlg.ShowModal() == wx.ID_OK:
            inv = dlg.inv.GetStringSelection()
            phi = float(dlg.angle.GetValue())
            offset = float(dlg.offset.GetValue())
            
            xinv, yinv = scroll_geo.coords_inv(phi, self.geo, 0, key_dict[inv])
            nxinv, nyinv = scroll_geo.coords_norm(phi, self.geo, 0, key_dict[inv])
            
            self.x.SetValue(str(xinv - nxinv[0]*offset))
            self.y.SetValue(str(yinv - nyinv[0]*offset))
            
            self.OnRefresh()
            
        dlg.Destroy()
        
    def OnRefresh(self, event = None):
        self.ax.cla()
        
        plotScrollSet(pi/4.0, 
                      axis = self.ax, 
                      geo = self.geo,
                      offsetScroll = self.geo.phi_ie_offset > 0)
        
        xlims = self.ax.get_xlim()
        ylims = self.ax.get_ylim()
        try:            
            x = float(self.x.GetValue())
            y = float(self.y.GetValue())
            if ylims[0] < y < ylims[1] and xlims[0] < x < xlims[1]:
                self.ax.plot(x,y,'yo')
        except ValueError:
            pass
        
        self.PP.canvas.draw()
        
class SuperButton(wx.Button):
    """ Button that destroys itself if right-clicked """
    def __init__(self, parent, *args, **kwargs):
        wx.Button.__init__(self, parent, *args, **kwargs)
        self.Bind(wx.EVT_RIGHT_UP, self.OnDestroy)
    def OnDestroy(self, event = None):
        self.Destroy()
        
class VirtualSensorsPanel(pdsim_panels.PDPanel):
    
    desc_map = dict()
    
    def __init__(self, parent, config, **kwargs):
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        # Now we are going to put everything into a scrolled window
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # The scrolled panel
        scrolled_panel = ScrolledPanel(self, size = (-1,-1), style = wx.TAB_TRAVERSAL, name="panel1")
        scrolled_panel.SetScrollbars(1, 1, 1, 1)
        
        self.AddSensor = wx.Button(scrolled_panel, label='Add Sensor')
        
        self.AddSensor.Bind(wx.EVT_BUTTON, self.OnAddSensor)
        
        self.sensor_sizer = wx.BoxSizer(wx.VERTICAL)
        
        s = textwrap.dedent('''
        INFORMATION: In this panel, virtual sensors can be added to
        the model.  These virtual sensors will "measure" all the state variables
        for the control volumes that overlap a given Cartesian coordinate.  In 
        this way it is possible to carry out virtual dynamic pressure measurements
        using the simulation code.  In this way you can check the selection of
        the locations of pressure sensors.
        
        You can add a sensor by clicking on the "Add Sensor" button below.  Sensors can
        be removed by right-clicking on the sensor in the list below''')
        
        self.description = wx.StaticText(scrolled_panel,label = s)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(HeaderStaticText(scrolled_panel, "Description"), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.description, 0, wx.ALIGN_CENTER_HORIZONTAL)

        sizer.Add(HeaderStaticText(scrolled_panel, "Virtual Sensors"), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.AddSensor, 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        sizer.Add(HeaderStaticText(scrolled_panel, "List of Virtual Sensors"), 0, wx.ALIGN_CENTER_HORIZONTAL)
        sizer.Add(self.sensor_sizer,0, wx.ALIGN_CENTER_HORIZONTAL)
        
        scrolled_panel.SetSizer(sizer)
        main_sizer.Add(scrolled_panel, 1, wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL)
        self.SetSizer(main_sizer)
        
        sizer.Layout()
        self.Fit()
        self.scrolled_panel = scrolled_panel
        
    def OnAddSensor(self, event):
        
        Scroll = self.Parent.panels_dict['GeometryPanel'].Scroll
        dlg = AddSensorDialog(None, Scroll.geo)
        if dlg.ShowModal() == wx.ID_OK:
            x,y = float(dlg.x.GetValue()), float(dlg.y.GetValue())
            but = SuperButton(self.scrolled_panel,label='x = {x:g}, y = {y:g}'.format(x=x,y=y))
            but.xval = x
            but.yval = y
            self.sensor_sizer.Add(but)
            self.sensor_sizer.Layout()
            self.GetSizer().Layout()
            self.Refresh()
            
        dlg.Destroy()
    
    def get_script_chunks(self, plugin_chunks = None):
        """ Chunk for the script file """
        chunk = ''
        for button in self.sensor_sizer.Children:
            x,y = button.Window.xval, button.Window.yval
            chunk += 'sim.add_sensor({x:g}, {y:g})\n'.format(x = x, y = y)
        return chunk
        
    def get_config_chunk(self):
        configdict = {}
        coords = []
        for button in self.sensor_sizer.Children:
            x,y = button.Window.xval, button.Window.yval
            coords.append((x,y))
        configdict['coords'] = coords
        return configdict