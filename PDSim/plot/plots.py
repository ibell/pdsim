from __future__ import print_function

import wx
import wx.aui
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as Toolbar
import numpy as np
from CoolProp.CoolProp import Props
import h5py

import matplotlib.pyplot as plt
from cycler import cycler

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
    def __init__(self,  Simulation, parent, id = -1, plot_names=None, family = None):
        wx.Panel.__init__(self, parent, id=id)
        self.nb = wx.aui.AuiNotebook(self)
        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Sim=Simulation
        self.family = family
        
        self.build_main_page()        
        
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
        page.figure.gca().set_prop_cycle(
            cycler('color', ['r', 'g', 'b', 'y', 'm', 'c']) *
            cycler('linestyle', ['-', '--', '-.'])
            )

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
                           ('Temperature-pressure',self.temperature_pressure),             
                           ('Heat transfer v. crank angle', self.heat_transfer),
                           ('Initial temperature history', self.initial_temperature_history),
                           ('Lump residuals v. lump temps', self.lumps_residual_v_lump_temps),
                           ('Discharge residual history', self.discharge_residual_history),
                           ]
        self.recip_plot_buttons = [('Valve lift v. crank angle',self.valve_theta)]
        self.scroll_plot_buttons = [('Pressure profile',self.pressure_profile),
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
        for value, callbackfcn in self.plot_buttons:
            btn = wx.Button(page, label = value)
            sizer.Add(btn)
            btn.Bind(wx.EVT_BUTTON, callbackfcn)
        
        if self.family == 'Scroll Compressor':
            more_plot_buttons = self.scroll_plot_buttons
        elif self.family == 'Recip Compressor':
            more_plot_buttons = self.recip_plot_buttons
        else:
            more_plot_buttons = None
        
        if more_plot_buttons is not None:           
            for value, callbackfcn in more_plot_buttons:
                btn = wx.Button(page, label = value)
                sizer.Add(btn)
                btn.Bind(wx.EVT_BUTTON, callbackfcn)
        else:
            print('could not add more buttons particular to current family:', self.family)
            
        page.SetSizer(sizer)
        self.nb.AddPage(page,"Main")
    
    def get_keys(self):
        if isinstance(self.Sim,h5py.File):
            NCV = self.Sim.get('/CVs/N').value
            keys = [self.Sim.get('/CVs/keys/'+str(i)).value for i in range(NCV)]
        else:
            keys = self.Sim.CVs.keys
        return keys

    def get(self, key):
        if isinstance(self.Sim, h5py.File):
            out = self.Sim.get('/' + key).value
        else:
            out=getattr(self.Sim, key)
        return out
        
    def stepsize_theta(self,event=None):
        #Stepsize
        axes = self.add('Stepsize').gca()
        theta = self.get('t')
        h = theta[1::]-theta[0:len(theta)-1]
        axes.semilogy(theta[0:len(theta)-1], h)
        axes.set_ylabel('Stepsize [rad]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def V_theta(self, event=None):
        # Volume
        axes = self.add('Volume').gca()
        theta, V = self.get('t'), self.get('V')
        V[V<1e-15]=np.nan
        axes.plot(theta, V.T, lw = 1.5)
        axes.set_ylabel('Volume [cm$^{3}$]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
    
    def dV_dtheta(self,event=None):
        #Derivative of Volume
        axes = self.add('Vol. Derivative').gca()
        theta, dV = self.get('t'), self.get('dV')
        dV *= 1e6
        dV[np.abs(dV)<1e-15]=np.nan
        axes.plot(theta, dV.T, lw = 1.5)
        axes.set_ylabel('Volume Derivative [cm$^{3}$/rad]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
    
    def T_theta(self,event=None):
        #Temperature
        axes = self.add('Temperature').gca()
        theta, T = self.get('t'), self.get('T')
        T[T<0.1]=np.nan
        axes.plot(theta, T.T, lw = 1.5)
        axes.set_ylabel('Temperature [K]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
        
    def p_theta(self,event=None):    
        #pressure
        axes = self.add('Pressure').gca()
        theta, p = self.get('t'), self.get('p')
        p[p<0.1] = np.nan
        axes.plot(theta, p.T, lw = 1.5)
        axes.set_ylabel('Pressure [kPa]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
        
    def p_V(self, event=None):
        #pressure-volume
        axes = self.add('P-V').gca()
        V, p = self.get('V'), self.get('p')
        V *= 1e6
        p[p<0.1] = np.nan
        V[V<1e-15] = np.nan
        axes.plot(V.T, p.T, lw = 1.5)
        axes.set_ylabel('Pressure [kPa]')
        axes.set_xlabel(r'Volume [cm$^{3}$]')
        
    def rho_theta(self,event=None):    
        #density
        axes = self.add('Density').gca()
        theta, rho = self.get('t'), self.get('rho')
        rho[rho<0.1]=np.nan
        axes.plot(theta, rho.T, lw = 1.5)
        axes.set_ylabel('Density [kg/m$^{3}$]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
        
    def m_theta(self,event=None):
        #Mass
        axes = self.add('Mass').gca()
        theta, rho, V = self.get('t'), self.get('rho'), self.get('V')
        m = rho*V
        m[m<1e-20]=np.nan
        axes.plot(theta, m.T, lw = 1.5)
        axes.set_ylabel('Mass [kg]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
        
    def mdot_theta(self,event=None):    
        #Mass Flow
        axes = self.add('Mdot').gca()
        axes.plot(self.Sim.FlowsProcessed.t,self.Sim.FlowsProcessed.summed_mdot['outlet.1'], lw = 1.5)
        axes.plot(self.Sim.FlowsProcessed.t,self.Sim.FlowsProcessed.summed_mdot['inlet.2'], lw = 1.5)
        axes.plot(self.Sim.FlowsProcessed.t,self.Sim.FlowsProcessed.mean_mdot['outlet.1']*np.ones_like(self.Sim.FlowsProcessed.t), lw = 1.5)
        axes.plot(self.Sim.FlowsProcessed.t,self.Sim.FlowsProcessed.mean_mdot['inlet.2']*np.ones_like(self.Sim.FlowsProcessed.t), lw = 1.5)
        axes.set_xlabel(r'$\theta$ [rad]')
        axes.set_ylabel(r'$\dot m$ [kg/s]')

    def initial_temperature_history(self, event = None):
        axes = self.add('Initial Temperature History').gca()
        if isinstance(self.Sim, h5py.File):
            TTT = self.Sim.get('/solvers/initial_states_history').value
        else:
            TTTT = self.Sim.solvers.initial_states_history
        xx = np.array(list(range(TTT.shape[0])))
        yy = TTT
        axes.plot(xx, yy, 'o-')
        axes.set_xlabel('Iteration Number')
        axes.set_ylabel('Temperature [K]')
        
    def lumps_residual_v_lump_temps(self, event = None):
        axes = self.add('Lump Error History').gca()
        if isinstance(self.Sim, h5py.File):
            Tlumps = self.Sim.get('/solvers/lump_eb_history/Tlumps').value
            lump_eb_error = self.Sim.get('/solvers/lump_eb_history/lump_eb_error').value
        else:
            Tlumps = self.solvers.lump_eb_history.Tlumps
            lump_eb_error = self.solvers.lump_eb_history.lump_eb_error
            
        axes.plot(Tlumps, lump_eb_error, 'o-')
        axes.set_xlabel('Lump temperature [K]')
        axes.set_ylabel('Lump energy balance [kW]')
        
    def discharge_residual_history(self, event = None):
        axes = self.add('Discharge Residual History').gca()
        if isinstance(self.Sim, h5py.File):
            hd = self.Sim.get('/solvers/hdisc_history/hd').value
            hd_error = self.Sim.get('/solvers/hdisc_history/hd_error').value
        else:
            hd = self.solvers.hdisc_history.hd
            hd_error = self.solvers.hdisc_history.hd_error
        axes.plot(hd, hd_error, 'o-')
        axes.set_xlabel('Discharge enthalpy [kJ/kg]')
        axes.set_ylabel('Discharge state error [kJ/kg]')
    
    def valve_theta(self, event = None):
        #valve lift
        if hasattr(self.Sim,'__hasValves__') and self.Sim.__hasValves__:
            axes = self.add('Valves').gca()
            axes.plot(self.Sim.t,self.Sim.xValves[0,:], lw = 1.5)
            axes.plot(self.Sim.t,self.Sim.xValves[2,:], lw = 1.5)
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
        Tsat = np.linspace(PropsSI(Fluid,'Tmin')+0.1, PropsSI(Fluid,'Tcrit')-1e-6)
        psat = np.array([PropsSI('P', 'T', T_, 'Q', 1.0, Fluid)/1000 for T_ in Tsat])
        axes.plot(Tsat, psat, lw = 1.5)
        axes.plot(T,p,'.')
        axes.set_xlabel('Temperature [K]')
        axes.set_ylabel(r'Pressure [kPa]')
            
    def heat_transfer(self, event = None):
        axes = self.add('Heat transfer').gca()
        theta, Q = self.get('t'), self.get('Q')
        Q[np.abs(Q)<1e-12]=np.nan
        axes.plot(theta, Q.T, lw = 1.5)
        axes.plot(theta,self.Sim.HTProcessed.summed_Q[0:self.Sim.Ntheta],lw=2)
        axes.set_ylabel(r'$\dot Q$ [kW]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
    
    def axial_force(self, event = None):
        #Axial force
            
        axes = self.add('Axial Force').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Fz = self.Sim.get('/forces/Fz').value.T
            summed_Fz = self.Sim.get('/forces/summed_Fz').value.T
        else:
            theta = self.Sim.t
            Fz = self.Sim.forces.Fz.T
            summed_Fz = self.Sim.forces.summed_Fz.T
        
        Fz[np.abs(Fz)<1e-12]=np.nan
        axes.plot(theta,Fz, lw = 1.5)
        axes.plot(theta,summed_Fz,lw=4)
        axes.set_ylabel(r'$F_z$ (only from the applied gas) [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
    
    def x_direction_force(self, event = None):
        #x-direction force
            
        axes = self.add('X Force').gca()
        
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Fx = self.Sim.get('/forces/Fx').value.T
        else:
            theta = self.Sim.t
            Fx = self.Sim.forces.Fx.T
            
        Fx[np.abs(Fx)<1e-12]=np.nan
        axes.plot(theta,Fx, lw = 1.5)
        axes.set_ylabel(r'$F_x$ [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
        
    def y_direction_force(self, event = None):
        #y-direction force

        axes = self.add('Y Force').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Fy = self.Sim.get('/forces/Fy').value.T
        else:
            theta = self.Sim.t
            Fy = self.Sim.forces.Fy.T
            
        Fy[np.abs(Fy)<1e-12]=np.nan
        axes.plot(theta,Fy, lw = 1.5)
        axes.set_ylabel(r'$F_y$ [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
    
    def force_component_trace(self, event = None):
        #trace of force components
            
        axes = self.add('Force trace').gca()
        if isinstance(self.Sim,h5py.File):
            Fx = self.Sim.get('/forces/Fx').value.T
            Fy = self.Sim.get('/forces/Fy').value.T
        else:
            Fx = self.Sim.forces.Fx.T
            Fy = self.Sim.forces.Fy.T
        axes.plot(Fx,Fy, lw = 1.5)
        axes.set_ylabel(r'$F_x$ [kN]')
        axes.set_xlabel(r'$F_y$ [kN]')            
        
    def force_trace(self, event = None):
        #trace of force components
            
        axes = self.add('Force trace').gca()

        if isinstance(self.Sim,h5py.File):
            Fx = self.Sim.get('/forces/summed_Fx').value.T
            Fy = self.Sim.get('/forces/summed_Fy').value.T
        else:
            Fx = self.Sim.forces.summed_Fx.T
            Fy = self.Sim.forces.summed_Fy.T
            
        Fx[np.abs(Fx)<1e-12]=np.nan
        Fy[np.abs(Fy)<1e-12]=np.nan
        
        axes.plot(Fx,Fy, lw = 1.5)
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
            Fm = self.Sim.forces.Fm.T
            
        Fm[np.abs(Fm)<1e-12] = np.nan
        axes.plot(theta, Fm, lw = 1.5)
        axes.set_ylabel(r'$F_m$ [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def radial_force(self, event = None):
        #Radial force magnitude
            
        axes = self.add('Radial force magnitude').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Fr = self.Sim.get('/forces/Fr').value.T
            mean_Fr = self.Sim.get('/forces/mean_Fr').value
            summed_Fr = self.Sim.get('/forces/summed_Fr').value
        else:
            theta = self.Sim.t
            Fr = self.Sim.forces.Fr.T
            mean_Fr = self.Sim.forces.mean_Fr
            summed_Fr = self.Sim.forces.summed_Fr
            
        Fr[np.abs(Fr)<1e-12]=np.nan
        axes.plot(theta,Fr, lw = 1.5)
        axes.plot(theta,mean_Fr*np.ones_like(theta),'k--', lw = 1.5)
        axes.plot(theta,summed_Fr,lw=4)
        axes.set_ylabel(r'$F_r$ (only from the applied gas) [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
        
    def tangential_force(self, event = None):
        #Tangential force magnitude
            
        axes = self.add('Tangential force magnitude').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            Ft = self.Sim.get('/forces/Ft').value.T
            mean_Ft = self.Sim.get('/forces/mean_Ft').value
            summed_Ft = self.Sim.get('/forces/summed_Ft').value
        else:
            theta = self.Sim.t
            Ft = self.Sim.forces.Ft.T
            mean_Ft = self.Sim.forces.mean_Ft
            summed_Ft = self.Sim.forces.summed_Ft
            
        Ft[np.abs(Ft)<1e-12]=np.nan
        axes.plot(theta,Ft, lw = 1.5)
        axes.plot(theta, mean_Ft*np.ones_like(theta), 'k--', lw = 1.5)
        axes.plot(theta,summed_Ft, lw = 4)
        axes.set_ylabel(r'$F_t$ (only from the applied gas) [kN]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(self.get_keys(), loc='upper right', bbox_to_anchor=(1, 1),
          ncol=2, fancybox=True, shadow=True)
    
        
    def torque(self, event = None):
        #Torque
        axes = self.add('Torque').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('/t').value
            tau = self.Sim.get('/forces/tau').value.T
            mean_tau = self.Sim.get('/forces/mean_tau').value
        else:
            theta = self.Sim.t
            tau = self.Sim.forces.tau.T
            mean_tau = self.Sim.forces.mean_tau
            
        tau[np.abs(tau)<1e-12]=np.nan
        axes.plot(theta,tau, lw = 1.5)
        axes.plot(theta, mean_tau*np.ones_like(theta),'k--', lw = 1.5)
        axes.set_ylabel(r'$\tau$ [kN-m]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
    def pressure_profile(self, event = None):
        #pressure profiles
            
        axes = self.add('Pressure profiles').gca()
        if isinstance(self.Sim,h5py.File):
            theta = self.Sim.get('summary/theta_profile').value
            p1 = self.Sim.get('summary/p1_profile').value.T
            p2 = self.Sim.get('summary/p2_profile').value.T
        else:
            theta = self.Sim.t
            p1 = self.Sim.summary.p1_profile
            p2 = self.Sim.summary.p2_profile
            
        axes.plot(theta, p1, lw = 1.5, label='s1-c1.x-d1-ddd')
        axes.plot(theta, p2, lw = 1.5, label='s2-c2.x-d2-ddd')
        axes.set_ylabel(r'Pressure [kPa]')
        axes.set_xlabel(r'$\theta$ [rad]')
        
        xmin,xmax = axes.get_xlim()
        axes.set_xlim(xmin, xmin + (xmax-xmin)*1.5)
        axes.legend(loc='upper right', bbox_to_anchor=(1, 1), ncol=2, 
                    fancybox=True, shadow=True)
         
    
def debug_plots(Comp, plotparent=None, plot_names = None):
    # Build a new frame, not embedded
    app = wx.App(False)
    frame = wx.Frame(None, -1,'Plotter')
    notebook = PlotNotebook(Comp, frame, plot_names=plot_names)
    frame.Show()
    app.MainLoop()

if __name__ == "__main__":
    print("File for doing plots.  Don't run this file")

