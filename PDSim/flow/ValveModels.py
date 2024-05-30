from __future__ import division, print_function

import numpy as np
from math import log,tan,sin,cos,pi

from PDSim.misc.datatypes import arraym, empty_arraym
from scipy.integrate import ode

# def public:
#     OUTPUT_VELOCITY
#     OUTPUT_MA

class ValveModel(object):
    """
    
    .. plot::

        import matplotlib.pyplot as plt
        import numpy as np
        
        fig=plt.figure()
        ax1=fig.add_subplot(111)
        ax1.fill(np.r_[0.2,0.4,0.4,0.2,0.2],np.r_[0,0,1,1,0],'grey')
        ax1.text(0.41,0.66,r'$\\leftarrow p_{low}A_{valve}$',size=20,ha='left',va='center')
        ax1.text(0.41,0.33,r'$\\leftarrow k_{valve}x$',size=20,ha='left',va='center')
        ax1.text(0.19,0.66,r'$p_{high}A_{valve}\\rightarrow$',size=20,ha='right',va='center')
        ax1.text(0.19,0.33,r'$\\frac{1}{2}C_D\\rho (V-\dot x)^2 A_{valve}\\rightarrow$',size=20,ha='right',va='center')
        ax1.set_xlim(0.1-1,0.1+1)
        ax1.axis('equal')
        ax1.axis('off')
        ax1.set_title('Pressure-dominant Free-body-diagram')
        plt.show()
    
    Pressure-dominant Regime
    
    .. math::
    
        M_{valve}\ddot x_{valve}+k_{valve}x_{valve} = (p_{high}-p_{low}) A_{valve}+\\frac{1}{2}C_D\\rho (V-\dot x_{valve})^2A_{valve}
        
    Two variables are :math:`x_2=\dot x_{valve}` and :math:`x_1=x_{valve}` where :math:`\ddot x_{valve}=\\frac{d}{dt}[\dot x_{valve}]` or :math:`\ddot x_{valve}=\dot x_2`.  Thus the system of derivatives is
    
    .. math::
    
        \mathbf f_{valves}=\\frac{d}{dt}\left[ \\begin{array}{c} \\dot x_{valve} \\\\ x_{valve} \end{array} \\right]=\left[ \\begin{array}{c} \\frac{d}{dt}[\dot x_{valve}] \\\\ \\frac{d}{dt}[x_{valve}] \end{array} \\right]
    
    Thus the system of equations is given by
    
    .. math::
    
        \dot x_1 = x_2
        
        M_{valve}\dot x_2+k_{valve}x_1 = (p_{high}-p_{low}) A_{valve}+\\frac{1}{2}C_D\\rho (V-x_2)^2A_{valve}
    
    which yields the solution for the derivatives of :math:`x_1` and :math:`x_2` of
    
    .. math::
    
        \dot x_1 = x_2
        
        \dot x_2 = \dfrac{(p_{high}-p_{low}) A_{valve}+\\frac{1}{2}C_D\\rho (V-x_2)^2A_{valve}-k_{valve}x_1}{M_{valve}}
        
    .. plot::
    
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig=plt.figure()
        ax1=fig.add_subplot(111)
        ax1.fill(np.r_[0.2,0.4,0.4,0.2,0.2],np.r_[0,0,1,1,0],'grey')
        ax1.text(0.41,0.66,r'$\leftarrow p_{low}A_{valve}$',size=20,ha='left',va='center')
        ax1.text(0.41,0.33,r'$\leftarrow k_{valve}x$',size=20,ha='left',va='center')
        ax1.text(0.19,0.8,r'$p_{low}A_{valve} \\rightarrow$',size=20,ha='right',va='center')
        ax1.text(0.19,0.5,r'$\\frac{1}{2}C_D\\rho (V-\dot x)^2 A_{valve} \\rightarrow$',size=20,ha='right',va='center')
        ax1.text(0.19,0.2,r'$\\rho (V-\dot x)^2 A_{port} \\rightarrow$',size=20,ha='right',va='center')
        ax1.set_xlim(0.1-1,0.1+1)
        ax1.axis('equal')
        ax1.axis('off')
        ax1.set_title('Mass-flux-dominant Free-body-diagram')
        plt.show()
    
    And if mass-flux-dominated, force balance given by
    
    .. math::
    
        M_{valve}\ddot x_{valve}+k_{valve}x_{valve} = \\frac{1}{2}C_D\\rho (\mathbf V-\dot x_{valve})^2 A_{valve}+\\rho (\mathbf V-\dot x_{valve})^2 A_{port}
    
    Which yields the solution for the system of derivatives of 
    
    .. math::
    
        \dot x_1 = x_2
        
        \dot x_2= \dfrac{\\frac{1}{2}C_D\\rho (\mathbf V-x_2)^2 A_{valve}+\\rho (\mathbf V-x_2)^2 A_{port}-k_{valve}x_1}{M_{valve}}
    """
    
    def __init__(self, d_valve, d_port, C_D, rho_valve, x_stopper,m_eff,k_valve,x_tr,key_up, key_down):
        '''
        # this function will be run when run the valve model first time
        '''
        self.xv = empty_arraym(2)
        self.rho_valve = rho_valve
        self.d_valve = d_valve
        self.d_port = d_port
        self.A_valve = pi*d_valve**2/4.0
        self.A_port = pi*d_port**2/4.0
        self.m_eff = m_eff
        self.C_D = C_D
        self.k_valve = k_valve
        self.x_stopper = x_stopper
        self.key_up = key_up
        self.key_down = key_down
        self.x_tr = x_tr 

    def get_States(self, Core):
        """
        Core is the main model core, it contains information that 
        is needed for the flow models
        
        self.CVs=ControlVolumeCollection()
        self.exists_keys = [self.keys[i] for i in self.exists_indices]
        self.exists_indices = [i for i in self.indices if self.CVs[i].exists]
        CVs : list of control volumes
        
        self.Tubes = TubeCollection()
        
        """
        exists_keys=Core.CVs.exists_keys    #        exists : bool ``True`` if control volume exists, ``False`` otherwise
        Tubes_Nodes=Core.Tubes.Nodes
              
        for key, Statevar in [(self.key_up,'State_up'),(self.key_down,'State_down')]:

            # Update the pointers to the states for the ends of the flow path
            if key in exists_keys:
                setattr(self,Statevar,Core.CVs[key].State)
            elif key in Tubes_Nodes:
                setattr(self,Statevar,Tubes_Nodes[key])  
        
    def set_xv(self,xv):

        self.xv = xv.copy()
        # if self.xv[0] < -1e-15 and self.xv.get_index(1) < 1e-15:
        if self.xv[0] < 0 and self.xv[1] < 1e-15:         
            self.xv[0] = 0.0
            self.xv[1] = 0.0
        #If it predicts a valve opening greater than max opening, just use the max opening
        elif self.xv[0] > self.x_stopper and self.xv[1] > 0: 
            self.xv[0] = self.x_stopper
            self.xv[1] = 0.0

    def _pressure_dominant(self, f, x,xdot,rho,V, deltap):
        """
        set_index() in cython becomes standard array assignment 
        """

        f[0] = xdot

        if abs(V-xdot) > 0:
            f[1] = ((V-xdot)/abs(V-xdot)*0.5*self.C_D*rho*V**2*self.A_valve+deltap*self.A_valve-self.k_valve*x)/(self.m_eff) #xdot
        else:
            f[1] = (deltap*self.A_valve-self.k_valve*x)/(self.m_eff) #d(xdot)dt

        return
        
    def _flux_dominant(self, f,x, xdot,rho,V):

        f[0] = xdot #dxdt

        if abs(V-xdot) > 0:
            f[1] = ((V-xdot)/abs(V-xdot)*0.5*self.C_D*rho*V**2*self.A_valve+(V-xdot)/abs(V-xdot)*rho*(V-xdot)**2*self.A_port-self.k_valve*x)/(self.m_eff)   #xdot
        else:
            f[1] = (-self.k_valve*x)/(self.m_eff) #d(xdot)dt
    
        return
        
    def A(self):

        if self.xv is None:
            print('self.xv is None')
        x = self.xv[0]
        if x >= self.x_tr:
            return self.A_port          # flux domain
        else:
            return pi*x*self.d_valve      # pressure domain
    
    def flow_velocity(self, State_up, State_down):   
        """
        For a given set of states, and a known valve lift, first
        check whether it is within the valve lift range, and then
        calculate the flow velocity
        """
        A = self.A()                # From def A to get the calculation area
        # print ' A in flow_velocity is:',A  # A is 0 all the time
        x = self.xv[0]
        
        if A > 0:
            if x > self.x_tr:
                return IsentropicNozzle(A, State_up, State_down, 2) # use the default mdot as our output
            else:
                # return x/self.x_tr*IsentropicNozzle(A, State_up, State_down, OUTPUT_VELOCITY)
                return x/self.x_tr*IsentropicNozzle(A, State_up, State_down, 2)    # can output velocity in cpython with public 
        else:
            return 0.0
        
    def derivs(self, Core):
        """
        Return the position and velocity as an arraym for the valve
        
        Parameters
        ----------
        Core : :class:`PDSimCore <PDSim.core.core.PDSimCore>` instance
        
        Returns
        -------
        out_array : :class:`arraym <PDSim.misc.datatypes.arraym>` instance
        
        """   
        f = empty_arraym(2)
        out_array = empty_arraym(2)
        x = self.xv[0]
        xdot = self.xv[1]
        
        self.get_States(Core)

        rho = self.State_up.get_rho()
        p_high = self.State_up.get_p()
        p_low = self.State_down.get_p()
        deltap = (p_high - p_low)*1000
        
        if deltap > 0:           
            V = self.flow_velocity(self.State_up, self.State_down)
        else:
            V = -self.flow_velocity(self.State_down, self.State_up)
        
        if x <= self.x_tr:
            self._pressure_dominant(f,x,xdot,rho,V,deltap)
        else:
            self._flux_dominant(f,x,xdot,rho,V)

        omega = Core.omega

        out_array[0] = f[0]/omega
        out_array[1] = f[1]/omega

        #TODO: check here (3)!!!
        if abs(x) < 1e-15 and xdot < -1e-12:
            # print 'stationary valve'
            out_array[0] = 0.0
            out_array[1] = 0.0
        
        return out_array #[dxdtheta, d(xdot)_dtheta]
    
    #TODO: check __cdict__
    def __dict__(self):
       
        # print 'Im in dict!'
        items=['A_port','A_valve','d_valve','h_valve','d_port','m_eff','C_D','a_valve','l_valve',
               'rho_valve','k_valve','x_stopper','key_up','key_down','xv','State_up', 'State_down',
               'E']
        
        return {item:getattr(self,item) for item in items}
    
    def __repr__(self):
        """
        Return a representation of Valve Model for outputting to screen
        """
        # print 'Im in repr!'
        s=''
        for item in self.__dict__():
            s += item+' : '+str(getattr(self,item))+'\n'
        return s
            
    def __reduce__(self):
        # print 'Im in reduce!'
        return rebuildValveModel,(self.__dict__(),) #TODO: check __cdict__


class ValveModel_kvalve(object):
    """
    Mathison M.M., PhD Thesis, Purdue University, 2011
    
    kvalve = a*exp(b*xv) + c
    
    """    
    def __init__(self, valve_type,d_valve, d_port, C_D,cheta, x_stopper,m_valve,x_tr,key_up, key_down):
        '''
        # this function will be run when run the valve model first time
        '''
        self.xv = empty_arraym(2)
        
        self.valve_type = valve_type
        self.d_valve = d_valve
        self.d_port = d_port
        self.A_valve = pi*d_valve**2/4.0
        self.A_port = pi*d_port**2/4.0
        self.m_valve = m_valve
        self.C_D = C_D
        self.cheta = cheta
        self.x_stopper = x_stopper
        self.key_up = key_up
        self.key_down = key_down
        self.x_tr = x_tr 
        
        self.Bv = pi*(self.d_port * 1.14)**2/ 4              #Effective area of the discharging port

    
    def get_States(self, Core):
        """
        Core is the main model core, it contains information that 
        is needed for the flow models
        
        self.CVs=ControlVolumeCollection()
        self.exists_keys = [self.keys[i] for i in self.exists_indices]
        self.exists_indices = [i for i in self.indices if self.CVs[i].exists]
        CVs : list of control volumes
        
        self.Tubes = TubeCollection()
        
        """
        exists_keys=Core.CVs.exists_keys    #        exists : bool ``True`` if control volume exists, ``False`` otherwise
        Tubes_Nodes=Core.Tubes.Nodes
              
        
        for key, Statevar in [(self.key_up,'State_up'),(self.key_down,'State_down')]:

            # Update the pointers to the states for the ends of the flow path
            if key in exists_keys:
                setattr(self,Statevar,Core.CVs[key].State)
            elif key in Tubes_Nodes:
                setattr(self,Statevar,Tubes_Nodes[key])  
        
    def set_xv(self,xv):

        self.xv = xv.copy()

        # if self.xv[0] < -1e-15 and self.xv.get_index(1) < 1e-15:
        if self.xv[0] <= 0:# and self.xv[1] < 1e-15:         
            self.xv[0] = 0.0
            self.xv[1] = 0.0
        #If it predicts a valve opening greater than max opening, just use the max opening
        elif self.xv[0] > self.x_stopper: # and self.xv[1] > 0: 
            self.xv[0] = self.x_stopper
            self.xv[1] = 0.0

    def _pressure_dominant(self, f, x,xdot,Dl,deltap):
        """
        set_index() in cython becomes standard array assignment 
        """
        f[0] = xdot     # Velovity of a valve
        f[1] = Dl - 2*self.cheta * self.omega_n * xdot - self.omega_n**2*x          #Acceleration of a valve

        return
        

    def A(self):

        if self.xv is None:
            print('self.xv is None')
        x = self.xv[0]
        if x <= 0.0:
            return 0.0
        elif x >= self.x_tr:
            return self.A_port          # flux domain
        else:
            Av = self.Bv / np.sqrt(1.5*(self.Bv / pi*x*self.d_valve)**2)            
        return Av       

    def derivs(self, Core):
        """
        Return the position and velocity as an arraym for the valve
        
        Parameters
        ----------
        Core : :class:`PDSimCore <PDSim.core.core.PDSimCore>` instance
        
        Returns
        -------
        out_array : :class:`arraym <PDSim.misc.datatypes.arraym>` instance
        
        """   
        f = empty_arraym(2)
        out_array = empty_arraym(2)
        x = self.xv[0]
        xdot = self.xv[1]        
        self.get_States(Core)
        rho = self.State_up.get_rho()
        p_high = self.State_up.get_p()
        p_low = self.State_down.get_p()
        deltap = (p_high - p_low)*1000
        self.m_vax = self.m_valve*(1 - 0.5*x/self.x_stopper)        #Mass of a valve moving as the free body
        
        if self.valve_type == '1stage':
            self.k_valve = 0.11584 *np.exp(x/0.00042) + 762.18456  
        
        elif self.valve_type == '2stage':
            self.k_valve = 0.11584*np.exp(x/0.00042) + 762.18456 
        
        self.F_fac = 0.04* (x/self.d_port-1)**2 + 0.06    #proportional factor based on the experiment
        self.omega_n = np.sqrt(self.k_valve/self.m_vax)   #Natural frequency
        Dl = self.F_fac*self.Bv*deltap/self.m_vax         #Acceleration of the gas againt the valve
        
        self._pressure_dominant(f,x,xdot,Dl,deltap)

        omega = Core.omega

        out_array[0] = f[0]/omega
        out_array[1] = f[1]/omega

        if abs(x) < 1e-15 and xdot < -1e-12:
            # print 'stationary valve'
            out_array[0] = 0.0
            out_array[1] = 0.0
        
        return out_array #[dxdtheta, d(xdot)_dtheta]
    
    def __dict__(self):
       
        # print 'Im in dict!'
        items=['A_port','A_valve','d_valve','d_port','m_valve','C_D',
               'k_valve','x_stopper','key_up','key_down','xv','State_up', 'State_down']
        
        return {item:getattr(self,item) for item in items}
    
    def __repr__(self):
        """
        Return a representation of Valve Model for outputting to screen
        """
        # print 'Im in repr!'
        s=''
        for item in self.__dict__():
            s += item+' : '+str(getattr(self,item))+'\n'
        return s
            
    def __reduce__(self):
        # print 'Im in reduce!'
        return rebuildValveModel,(self.__dict__(),) 


########################################################
# Add isentropicNozzle model into valve model
########################################################

def  IsentropicNozzle(A, State_up, State_down, other_output = -1):
    """
    The mass flow rate is calculated by using isentropic flow model
    
    Parameters
    ----------
    
    A : double
        Throat area of the nozzle [m\ :math:`^2`\ ]
    State_up : :class:`State <CoolProp.State.State>` instance
        Upstream ``State`` instance
    State_down : :class:`State <CoolProp.State.State>` instance
        Downstream ``State`` instance
    other_output : int
        Default is to return the mass flow rate, can over-ride by passing ``flow_models.OUTPUT_VELOCITY`` or ``flow_models.OUTPUT_MA`` instead
        
    Returns
    -------
    out : double
        Default is to return the mass flow rate, can over-ride by passing flags in the other_output variable
    
    """
    
    #some cython type declarations
    # cython.declare(T_up = cython.double, T_down = cython.double, mdot = cython.double)
    # Since ideal, R=cp-cv, and k=cp/cv
    cp = State_up.get_cp0()
    R = 8314.472/State_up.get_MM() #[J/kg/K]
    cv = cp-R/1000.0 #[kJ/kg/K]
    k = cp / cv
    
    p_up = State_up.get_p()
    T_up = State_up.get_T()
    p_down = State_down.get_p()
    
    # Speed of sound
    c = (k*R*T_up)**0.5
    # Upstream density
    rho_up = p_up*1000.0/(R*T_up)
    pr = p_down/p_up
    pr_crit = (1+(k-1)/2)**(k/(1-k))
    
    if pr > pr_crit: 
        # Mass flow rate if not choked [kg/s]
        mdot = A*p_up*1000.0/(R*T_up)**0.5*(2*k/(k-1.0)*pr**(2.0/k)*(1-pr**((k-1.0)/k)))**0.5
        # Throat temperature [K]
        T_down = T_up*(p_down/p_up)**((k-1.0)/k)
        # Throat density [kg/m3]
        rho_down = p_down*1000.0/(R*T_down)
        # Velocity at throat
        v = mdot/(rho_down*A)
        # Mach number
        Ma = v/c
    else:
        # Mass flow rate if choked
        mdot = A*rho_up*(k*R*T_up)**0.5*(1.+(k-1.)/2.)**((1+k)/(2*(1-k)))
        # Velocity at throat
        v = c
        # Mach Number
        Ma = 1.0
    
    if other_output < 0:
        return mdot
    # elif other_output == OUTPUT_VELOCITY:
    elif other_output == 2:
        return v
    elif other_output == OUTPUT_MA:
        return Ma
   

    
    