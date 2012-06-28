
import wx
import wx.aui
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar
import numpy as np
from CoolProp.CoolProp import Props

class Plot(wx.Panel):
    def __init__(self, parent, id = -1, dpi = None, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)
        self.figure = mpl.figure.Figure(dpi=dpi, figsize=(2,2))
        self.canvas = Canvas(self, -1, self.figure)
        self.toolbar = Toolbar(self.canvas)
        self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)
        sizer.Add(self.toolbar, 0 , wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)

class PlotNotebook(wx.Panel):
    def __init__(self,  Simulation, parent, id = -1):
        wx.Panel.__init__(self, parent, id=id)
        self.nb = wx.aui.AuiNotebook(self)
        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.build_main_page()
        self.Sim=Simulation

    def update(self,Simulation):
        """
        Pass in a new Simulation instance
        
        Will clear all plots
        """
        #Update the structure
        self.Sim = Simulation
        while(self.nb.GetPageCount()>1):
            self.nb.DeletePage(1)
        
    def add(self,name="plot"):
        page = Plot(self.nb)
        self.nb.AddPage(page,name)
        return page.figure
    
    def build_main_page(self):
        page = wx.Panel(self.nb,-1)
        sizer = wx.BoxSizer(wx.VERTICAL)
        label1=wx.StaticText(page,label='Click on the buttons below to add plot')
        sizer.Add(label1)
        plot_buttons=[('Volume v. crank angle',self.V_theta),
                      ('Derivative of Volume v. crank angle',self.dV_dtheta),
                      ('Temperature v. crank angle',self.T_theta),
                      ('Pressure v. crank angle',self.p_theta),
                      ('Pressure v. volume',self.p_V),
                      ('Density v. crank angle',self.rho_theta),
                      ('Mass v. crank angle',self.m_theta),
                      ('Mass flow v. crank angle',self.mdot_theta),
                      ('Valve lift v. crank angle',self.valve_theta)
                      ]
        for value,callbackfcn in plot_buttons:
            btn = wx.Button(page,label=value)
            sizer.Add(btn)
            btn.Bind(wx.EVT_BUTTON,callbackfcn)
        page.SetSizer(sizer)
        self.nb.AddPage(page,"Main")

    def V_theta(self,event=None):
        #Volume
        axes = self.add('Volume').gca()
        theta=self.Sim.t
        V=self.Sim.V.T*1e6
        V[V<1e-15]=np.nan
        axes.plot(theta,V)
        axes.set_ylabel('Volume [cm$^{3}$]')
        axes.set_xlabel(r'$\theta$ [rad]')
    
    def dV_dtheta(self,event=None):
        #Derivative of Volume
        axes = self.add('Vol. Derivative').gca()
        theta=self.Sim.t
        dV=self.Sim.dV.T*1e6
        dV[np.abs(dV)<1e-15]=np.nan
        axes.plot(theta,dV)
        axes.set_ylabel('Volume Derivative [cm$^{3}$/rad]')
        axes.set_xlabel(r'$\theta$ [rad]') 
    
    def T_theta(self,event=None):
        #Temperature
        axes = self.add('Temperature').gca()
        theta=self.Sim.t
        T=self.Sim.T.T
        T[T<0.1]=np.nan
        axes.plot(theta,T)
        axes.set_ylabel('Temperature [K]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def p_theta(self,event=None):    
        #pressure
        axes = self.add('Pressure').gca()
        theta=self.Sim.t
        p=self.Sim.p.T
        p[p<0.1]=np.nan
        axes.plot(theta,p)
        axes.set_ylabel('Pressure [kPa]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def p_V(self, event=None):
        #pressure-volume
        axes = self.add('P-V').gca()
        p=self.Sim.p.T
        p[p<0.1]=np.nan
        V=self.Sim.V.T*1e6
        V[V<1e-15]=np.nan
        axes.plot(V,p)
        axes.set_ylabel('Pressure [kPa]')
        axes.set_xlabel(r'Volume [cm$^{3}$]')
        
    def rho_theta(self,event=None):    
        #density
        axes = self.add('Density').gca()
        theta=self.Sim.t
        rho=self.Sim.rho.T
        rho[rho<0.1]=np.nan
        axes.plot(theta,rho)
        axes.set_ylabel('Density [kg/m$^{3}$]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def m_theta(self,event=None):
        #Mass
        axes = self.add('Mass').gca()
        theta=self.Sim.t
        m=self.Sim.rho.T*self.Sim.V.T
        m[m<1e-20]=np.nan
        axes.plot(theta,m)
        axes.set_ylabel('Mass [kg]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def mdot_theta(self,event=None):    
        #Mass Flow
        axes = self.add('Mdot').gca()
        axes.plot(self.Sim.FlowsProcessed.t,self.Sim.FlowsProcessed.summed_mdot['outlet.1'])
        axes.plot(self.Sim.FlowsProcessed.t,self.Sim.FlowsProcessed.summed_mdot['inlet.2'])
        axes.plot(self.Sim.FlowsProcessed.t,self.Sim.FlowsProcessed.mean_mdot['outlet.1']*np.ones_like(self.Sim.FlowsProcessed.t))
        axes.plot(self.Sim.FlowsProcessed.t,self.Sim.FlowsProcessed.mean_mdot['inlet.2']*np.ones_like(self.Sim.FlowsProcessed.t))
        axes.set_xlabel(r'$\theta$ [rad]')
        axes.set_ylabel(r'$\dot m$ [kg/s]')
    
    def valve_theta(self, event = None):
        #valve lift
        if hasattr(self.Sim,'__hasValves__') and self.Sim.__hasValves__:
            axes = self.add('Valves').gca()
            print self.Sim.xValves.shape
            axes.plot(self.Sim.t,self.Sim.xValves[0,:])
            axes.plot(self.Sim.t,self.Sim.xValves[2,:])
            axes.set_xlabel(r'$\theta$ [rad]')
            axes.set_ylabel(r'Valve lift [m]')
        
    #    if not Comp.__hasLiquid__:
    #        #Fluid T-p plot
    #        axes = notebook.add('T-P phase').gca()
    #        Fluid=Comp.CVs['s1'].State.Fluid
    #        #Saturation curve
    #        Tsat=np.linspace(Props(Fluid,'Ttriple'),Props(Fluid,'Tcrit'))
    #        psat=np.array([Props('P','T',T_,'Q',1.0,Fluid) for T_ in Tsat])
    #        axes.plot(Tsat,psat)
    #        axes.plot(T,p,'.')
    #        axes.set_xlabel('Temperature [K]')
    #        axes.set_ylabel(r'Pressure [kPa]')
    
def debug_plots(Comp,plotparent=None):
    #Build a new frame, not embedded
    app = wx.PySimpleApp()
    frame = wx.Frame(None,-1,'Plotter')
    notebook = PlotNotebook(Comp,frame)
    frame.Show()
    app.MainLoop()

if __name__ == "__main__":
    print "File for doing plots.  Don't run this file"

