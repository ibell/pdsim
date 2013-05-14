
import wx
import wx.aui
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar
import numpy as np
from CoolProp.CoolProp import Props
import h5py

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
    def __init__(self,  Simulation, parent, id = -1,plot_names=None):
        wx.Panel.__init__(self, parent, id=id)
        self.nb = wx.aui.AuiNotebook(self)
        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.build_main_page()
        self.Sim=Simulation
        
        if plot_names is not None:
            for plot_name in plot_names:
                matched = False
                for name, func in self.plot_buttons:
                    if name == plot_name:
                        func()
                        matched = True 
                #Not matched, quit
                if not matched:
                    raise ValueError("Could not match the name: '"+plot_name+'"')

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
        sizer = wx.FlexGridSizer(cols=2)
        label1=wx.StaticText(page,label='Click on the buttons below to add plot')
        sizer.Add(label1)
        self.plot_buttons=[('Stepsize',self.stepsize_theta),
                           ('Volume v. crank angle',self.V_theta),
                          ('Derivative of Volume v. crank angle',self.dV_dtheta),
                          ('Temperature v. crank angle',self.T_theta),
                          ('Pressure v. crank angle',self.p_theta),
                          ('Pressure v. volume',self.p_V),
                          ('Density v. crank angle',self.rho_theta),
                          ('Mass v. crank angle',self.m_theta),
                          ('Mass flow v. crank angle',self.mdot_theta),
                          ('Valve lift v. crank angle',self.valve_theta),
                          ('Temperature-pressure',self.temperature_pressure),
                          ('Heat transfer v. crank angle', self.heat_transfer),
                          ('Axial force v. crank angle',self.axial_force),
                          ('X-direction force v. crank angle',self.x_direction_force),
                          ('Y-direction force v. crank angle',self.y_direction_force),
                          ('Crank pin force magnitude v. crank angle',self.magnitude_force),
                          ('Gas Torque v. crank angle',self.torque),
                          ('Force trace', self.force_trace),
                          ('Force component trace',self.force_component_trace),
                          ('Radial force', self.radial_force),
                          ('Tangential force', self.tangential_force)
                          
                          ]
        for value,callbackfcn in self.plot_buttons:
            btn = wx.Button(page,label=value)
            sizer.Add(btn)
            btn.Bind(wx.EVT_BUTTON,callbackfcn)
        page.SetSizer(sizer)
        self.nb.AddPage(page,"Main")
    
    def stepsize_theta(self,event=None):
        #Stepsize
        axes = self.add('Stepsize').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
        else:
            theta=self.Sim.t
        h = theta[1::]-theta[0:len(theta)-1]
        axes.semilogy(theta[0:len(theta)-1],h)
        axes.set_ylabel('Stepsize [rad]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def V_theta(self,event=None):
        #Volume
        axes = self.add('Volume').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            V = self.Sim.get('/V').value.T*1e6
        else:
            theta=self.Sim.t
            V=self.Sim.V.T*1e6
        V[V<1e-15]=np.nan
        axes.plot(theta,V)
        axes.set_ylabel('Volume [cm$^{3}$]')
        axes.set_xlabel(r'$\theta$ [rad]')
    
    def dV_dtheta(self,event=None):
        #Derivative of Volume
        axes = self.add('Vol. Derivative').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            dV = self.Sim.get('/dV').value.T*1e6
        else:
            theta=self.Sim.t
            dV=self.Sim.dV.T*1e6
        dV[np.abs(dV)<1e-15]=np.nan
        axes.plot(theta,dV)
        axes.set_ylabel('Volume Derivative [cm$^{3}$/rad]')
        axes.set_xlabel(r'$\theta$ [rad]')
    
    def T_theta(self,event=None):
        #Temperature
        axes = self.add('Temperature').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            T = self.Sim.get('/T').value.T
        else:
            theta = self.Sim.t
            T = self.Sim.T.T
        T[T<0.1]=np.nan
        axes.plot(theta,T)
        axes.set_ylabel('Temperature [K]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def p_theta(self,event=None):    
        #pressure
        axes = self.add('Pressure').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            p = self.Sim.get('/p').value.T
        else:
            theta = self.Sim.t
            p = self.Sim.p.T
        p[p<0.1]=np.nan
        axes.plot(theta,p)
        axes.set_ylabel('Pressure [kPa]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def p_V(self, event=None):
        #pressure-volume
        axes = self.add('P-V').gca()
        if isinstance(self.Sim,h5py.File):
            V = self.Sim.get('/V').value.T*1e6
            p = self.Sim.get('/p').value.T
        else:
            V = self.Sim.V.T*1e6
            p = self.Sim.p.T
        p[p<0.1]=np.nan
        V[V<1e-15]=np.nan
        axes.plot(V,p)
        axes.set_ylabel('Pressure [kPa]')
        axes.set_xlabel(r'Volume [cm$^{3}$]')
        
    def rho_theta(self,event=None):    
        #density
        axes = self.add('Density').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            rho = self.Sim.get('/rho').value.T
        else:
            theta = self.Sim.t
            rho = self.Sim.rho.T
            
        rho[rho<0.1]=np.nan
        axes.plot(theta,rho)
        axes.set_ylabel('Density [kg/m$^{3}$]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def m_theta(self,event=None):
        #Mass
        axes = self.add('Mass').gca()
        
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            m = self.Sim.get('/rho').value.T*self.Sim.get('/V').value.T
        else:
            theta = self.Sim.t
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
            axes.plot(self.Sim.t,self.Sim.xValves[0,:])
            axes.plot(self.Sim.t,self.Sim.xValves[2,:])
            axes.set_xlabel(r'$\theta$ [rad]')
            axes.set_ylabel(r'Valve lift [m]')
        
    def temperature_pressure(self, event = None):
        #Fluid T-p plot
        axes = self.add('T-P phase').gca()
        if isinstance(self.Sim,h5py.File):
            Fluid = self.Sim.get('/Tubes/0/State1/Fluid').value
            T = self.Sim.get('/T').value.T
            p = self.Sim.get('/p').value.T
        else:
            Fluid = self.Sim.CVs.exists_CV[0].State.Fluid
            T = self.Sim.T
            p = self.Sim.p
            
        #Saturation curve
        Tsat = np.linspace(Props(Fluid,'Tmin')+0.1, Props(Fluid,'Tcrit')-1e-6)
        psat = np.array([Props('P', 'T', T_, 'Q', 1.0, Fluid) for T_ in Tsat])
        axes.plot(Tsat, psat)
        axes.plot(T,p,'.')
        axes.set_xlabel('Temperature [K]')
        axes.set_ylabel(r'Pressure [kPa]')
            
    def heat_transfer(self, event = None):
        #Axial force
            
        axes = self.add('Heat transfer').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Q = self.Sim.get('/Q').value.T
        else:
            theta = self.Sim.t
            Q = self.Sim.Q.T
            
        Q[np.abs(Q)<1e-12]=np.nan
        axes.plot(theta,Q)
        axes.plot(theta,self.Sim.HTProcessed.summed_Q[0:self.Sim.Ntheta],lw=2)
        axes.set_ylabel(r'$\dot Q$ [kW]')
        axes.set_xlabel(r'$\theta$ [rad]')
    
    def axial_force(self, event = None):
        #Axial force
            
        axes = self.add('Axial Force').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Fz = self.Sim.get('/forces/Fz').value.T
        else:
            theta = self.Sim.t
            Fz = self.Sim.Fz.T
        
        Fz[np.abs(Fz)<1e-12]=np.nan
        axes.plot(theta,Fz)
        axes.set_ylabel(r'$F_z$ [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
    
    def x_direction_force(self, event = None):
        #x-direction force
            
        axes = self.add('X Force').gca()
        
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Fx = self.Sim.get('/forces/Fx').value.T
        else:
            theta = self.Sim.t
            Fx = self.Sim.Fx.T
            
        Fx[np.abs(Fx)<1e-12]=np.nan
        axes.plot(theta,Fx)
        axes.set_ylabel(r'$F_x$ [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def y_direction_force(self, event = None):
        #y-direction force

        axes = self.add('Y Force').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Fy = self.Sim.get('/forces/Fy').value.T
        else:
            theta = self.Sim.t
            Fy = self.Sim.Fy.T
            
        Fy[np.abs(Fy)<1e-12]=np.nan
        axes.plot(theta,Fy)
        axes.set_ylabel(r'$F_y$ [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
    
    def force_component_trace(self, event = None):
        #trace of force components
            
        axes = self.add('Force trace').gca()
        if isinstance(self.Sim,h5py.File):
            Fx = self.Sim.get('/forces/Fx').value.T
            Fy = self.Sim.get('/forces/Fy').value.T
        else:
            Fx = self.Sim.Fx.T
            Fy = self.Sim.Fy.T
        axes.plot(Fx,Fy)
        axes.set_ylabel(r'$F_x$ [kN]')
        axes.set_xlabel(r'$F_y$ [kN]')    
        
    def force_trace(self, event = None):
        #trace of force components
            
        axes = self.add('Force trace').gca()

        if isinstance(self.Sim,h5py.File):
            Fx = self.Sim.get('/forces/summed_Fx').value.T
            Fy = self.Sim.get('/forces/summed_Fy').value.T
        else:
            Fx = self.Sim.summed_Fx.T
            Fy = self.Sim.summed_Fy.T
            
        Fx[np.abs(Fx)<1e-12]=np.nan
        Fy[np.abs(Fy)<1e-12]=np.nan
        
        axes.plot(Fx,Fy)
        axes.set_ylabel(r'$F_x$ [kN]')
        axes.set_xlabel(r'$F_y$ [kN]')
        
    def magnitude_force(self, event = None):
        #Crank pin force magnitude
            
        axes = self.add('Shaft force magnitude').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Fm = self.Sim.get('/forces/Fm').value.T
        else:
            theta = self.Sim.t
            Fm = self.Sim.Fm.T
            
        Fm[np.abs(Fm)<1e-12] = np.nan
        axes.plot(theta, Fm)
        axes.set_ylabel(r'$F_m$ [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def radial_force(self, event = None):
        #Radial force magnitude
            
        axes = self.add('Radial force magnitude').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Fr = self.Sim.get('/forces/Fr').value.T
            mean_Fr = self.Sim.get('/forces/mean_Fr').value
        else:
            theta = self.Sim.t
            Fr = self.Sim.Fr.T
            mean_Fr = self.Sim.forces.mean_Fr
            
        Fr[np.abs(Fr)<1e-12]=np.nan
        axes.plot(theta,Fr)
        axes.plot(theta,mean_Fr*np.ones_like(theta),'k--')
        axes.set_ylabel(r'$F_r$ [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def tangential_force(self, event = None):
        #Tangential force magnitude
            
        axes = self.add('Tangential force magnitude').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Ft = self.Sim.get('/forces/Ft').value.T
            mean_Ft = self.Sim.get('/forces/mean_Ft').value
        else:
            theta = self.Sim.t
            Ft = self.Sim.Ft.T
            mean_Ft = self.Sim.forces.mean_Ft
            
        Ft[np.abs(Ft)<1e-12]=np.nan
        axes.plot(theta,Ft)
        axes.plot(theta, mean_Ft*np.ones_like(theta), 'k--')
        axes.set_ylabel(r'$F_t$ [kN]')
        axes.set_xlabel(r'$\theta$ [rad]') 
    
        
    def torque(self, event = None):
        #Torque
        axes = self.add('Torque').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            tau = self.Sim.get('/forces/tau').value.T
            mean_tau = self.Sim.get('/forces/mean_tau').value
        else:
            theta = self.Sim.t
            tau = self.Sim.tau.T
            mean_tau = self.Sim.forces.mean_tau
            
        tau[np.abs(tau)<1e-12]=np.nan
        axes.plot(theta,tau)
        axes.plot(theta, mean_tau*np.ones_like(theta),'k--')
        axes.set_ylabel(r'$\tau$ [kN-m]')
        axes.set_xlabel(r'$\theta$ [rad]')
         
    
def debug_plots(Comp, plotparent=None, plot_names = None):
    #Build a new frame, not embedded
    app = wx.PySimpleApp()
    frame = wx.Frame(None,-1,'Plotter')
    notebook = PlotNotebook(Comp,frame,plot_names=plot_names)
    frame.Show()
    app.MainLoop()

if __name__ == "__main__":
    print "File for doing plots.  Don't run this file"

