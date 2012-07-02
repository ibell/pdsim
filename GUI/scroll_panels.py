# -*- coding: latin-1 -*-
import pdsim_panels
import wx
import numpy as np
from math import pi
from PDSim.scroll.core import Scroll
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as WXCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as WXToolbar
from PDSim.scroll.plots import plotScrollSet

LabeledItem = pdsim_panels.LabeledItem

class PlotPanel(wx.Panel):
    def __init__(self, parent, **kwargs):
        wx.Panel.__init__(self, parent, size = (300,200), **kwargs)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.figure = mpl.figure.Figure(dpi=100, figsize=(2, 2))
#        self.figure.set_figwidth(2.0)
#        self.figure.set_figheight(2.0)
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
    The geometry panel of the reciprocating compressor
    Loads all parameters from the configuration file
    """
    def __init__(self,parent,configfile,**kwargs):
        pdsim_panels.PDPanel.__init__(self,parent,**kwargs)
        
        #Loads all the parameters from the config file
        configdict, descdict = self.get_from_configfile('GeometryPanel')
        
        # Things in self.items are linked through to the module code where 
        # it attempts to set the attribute.  They are also automatically
        # written to configuration file
        self.items = [
        dict(attr='Vdisp'),
        dict(attr='Vratio'),
        dict(attr='t'),
        dict(attr='ro')
        ]
        
        sizerInputs = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        
        self.ConstructItems(self.items, sizerInputs, configdict, descdict)
        
        for item in self.items:
            setattr(self,item['attr'],item['textbox'])
        
        kwargs = dict(label = "phi_i0 [radian]",
                      tooltip = 'Initial involute angle for inner involute'
                      )
        self.phi_i0_label, self.phi_i0 = LabeledItem(self, **kwargs)
        
        self.phi_is_label, self.phi_is= LabeledItem(self,
                                                       label="phi_is [radian]")
        self.phi_ie_label, self.phi_ie= LabeledItem(self,
                                                       label="phi_ie [radian]")
        self.phi_o0_label, self.phi_o0= LabeledItem(self,
                                                       label="phi_o0 [radian]")
        self.phi_os_label, self.phi_os= LabeledItem(self,
                                                       label="phi_os [radian]")
        self.phi_oe_label, self.phi_oe= LabeledItem(self,
                                                       label="phi_oe [radian]")
        self.rb_label, self.rb = LabeledItem(self,
                                             label="rb [m]")
        self.hs_label, self.hs= LabeledItem(self,
                                            label="hs [m]")
        
        self.phi_i0.Enable(False)
        self.phi_is.Enable(False)
        self.phi_ie.Enable(False)
        self.phi_o0.Enable(False)
        self.phi_os.Enable(False)
        self.phi_oe.Enable(False)
        self.rb.Enable(False)
        self.hs.Enable(False)
        
        self.PP = PlotPanel(self)
        
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizerInputs)
        sizer.AddSpacer(10)
        sizer.Add(self.PP,0,wx.EXPAND)
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
        
        self.SetSizer(sizer)
        sizer.Layout()
        
        for item in self.items:
            item['textbox'].Bind(wx.EVT_KILL_FOCUS,self.OnChangeParam)
        
        # Keep a local copy of the scroll in order to be able to use the 
        # set_scroll_geo and set_disc_geo functions
        self.Scroll=Scroll()
        
        self.OnChangeParam()
        
    def skip_list(self):
        """
        Returns a list of atttributes to skip setting in calculate() function
        """
        return ['Vdisp','Vratio','t','ro']
        
    def OnChangeParam(self, event = None):
        
        Vdisp=float(self.Vdisp.GetValue())
        Vratio=float(self.Vratio.GetValue())
        t=float(self.t.GetValue())
        ro=float(self.ro.GetValue())
        
        self.Scroll.set_scroll_geo(Vdisp,Vratio,t,ro) #Set the scroll wrap geometry
        self.Scroll.set_disc_geo('2Arc', r2='PMP')
        
        self.phi_i0.SetValue(str(self.Scroll.geo.phi_i0))
        self.phi_is.SetValue(str(self.Scroll.geo.phi_is))
        self.phi_ie.SetValue(str(self.Scroll.geo.phi_ie))
        self.phi_o0.SetValue(str(self.Scroll.geo.phi_o0))
        self.phi_os.SetValue(str(self.Scroll.geo.phi_os))
        self.phi_oe.SetValue(str(self.Scroll.geo.phi_oe))
        self.rb.SetValue(str(self.Scroll.geo.rb))
        self.hs.SetValue(str(self.Scroll.geo.h))
        
        self.ax.cla()
        plotScrollSet(pi/4.0, axis = self.ax, geo = self.Scroll.geo)
        self.PP.canvas.draw()
        
class MassFlowPanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, configfile,**kwargs):
    
        pdsim_panels.PDPanel.__init__(self, parent,**kwargs)
        
        #Loads all the parameters from the config file
        self.configdict,self.descdict = self.get_from_configfile('MassFlowPanel')
        
        self.items1 = [
        dict(attr='d_discharge'),
        dict(attr='d_suction'),
        ]
        box_sizer = wx.BoxSizer(wx.VERTICAL)
        box_sizer.Add(wx.StaticText(self,-1,"Required Inputs"))
        box_sizer.Add(wx.StaticLine(self,-1,(25, 50), (300,1)))
        
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        self.ConstructItems(self.items1,sizer,self.configdict,self.descdict)

        box_sizer.Add(sizer)
        box_sizer.Add((10,10))
        box_sizer.Add(wx.StaticText(self,-1,"Valve Inputs"))
        box_sizer.Add(wx.StaticLine(self,-1,(25, 50), (300,1)))
        
        self.items2 = [
        dict(attr='valve_E'),
        dict(attr='valve_d'),
        dict(attr='valve_h'),
        dict(attr='valve_l'),
        dict(attr='valve_a'),
        dict(attr='valve_x_stopper'),
        dict(attr='valve_rho'),
        dict(attr='valve_C_D'),
        ]
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        self.ConstructItems(self.items2,sizer,self.configdict,self.descdict)
        box_sizer.Add(sizer)
        self.SetSizer(box_sizer)
        box_sizer.Layout()
        
        self.items=self.items1+self.items2

        
class MechanicalLossesPanel(pdsim_panels.PDPanel):
    
    def __init__(self, parent, configfile,**kwargs):
    
        pdsim_panels.PDPanel.__init__(self, parent, **kwargs)
        
        #Loads all the parameters from the config file (case-sensitive)
        self.configdict, self.descdict = self.get_from_configfile('MechanicalLossesPanel')
        
        self.items = [
        dict(attr='eta_motor'),
        dict(attr='h_shell'),
        dict(attr='A_shell'),
        dict(attr='Tamb'),
        dict(attr='mu_oil'),
        dict(attr='delta_gap'),
        ]
        
        sizer = wx.FlexGridSizer(cols=2, vgap=4, hgap=4)
        
        self.ConstructItems(self.items,sizer,self.configdict,self.descdict)

        self.SetSizer(sizer)
        sizer.Layout()