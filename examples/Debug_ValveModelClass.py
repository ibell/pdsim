from __future__ import division

import numpy as np
from math import log,tan,sin,cos,pi

from PDSim.misc.datatypes import arraym, empty_arraym




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
    
#     #Give Cython type definitions for terms
#     @cython.locals(d_valve=cython.double,d_port=cython.double,C_D=cython.double,
#                    h_valve=cython.double,a_valve=cython.double, l_valve=cython.double, 
#                    rho_valve=cython.double,E=cython.double,x_stopper=cython.double,
#                    key_up=cython.bytes,key_down=cython.bytes)
    def __init__(self, d_valve, d_port, C_D, h_valve, a_valve,
                 l_valve, rho_valve, E, x_stopper, key_up, key_down):
        
        
        I=(d_valve*h_valve**3)/12  #Moment of Intertia for discharge valve,[m^4]
        k_valve=(6*E*I)/(a_valve**2*(3*l_valve-a_valve))    #Valve stiffness
        m_eff=(1/3)*rho_valve*l_valve*d_valve*h_valve      #Effective mass of valve reeds
        
        self.E = E
        self.rho_valve = rho_valve
        self.a_valve = a_valve
        self.l_valve = l_valve
        self.h_valve = h_valve
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
        self.x_tr = 0.25*(self.d_port**2/self.d_valve)
        
        self.xv = empty_arraym(2)
        
    def get_States(self, Core):
        """
        Core is the main model core, it contains information that 
        is needed for the flow models
        """
        exists_keys=Core.CVs.exists_keys
        Tubes_Nodes=Core.Tubes.Nodes
        for key, Statevar in [(self.key_up,'State_up'),(self.key_down,'State_down')]:
            ## Update the pointers to the states for the ends of the flow path
            if key in exists_keys:
                setattr(self,Statevar,Core.CVs[key].State)
            elif key in Tubes_Nodes:
                setattr(self,Statevar,Tubes_Nodes[key])              
    
    #@cython.cdivision(True)
    def _pressure_dominant(self, f, x,xdot,rho,V, deltap):
        """
        set_index() in cython becomes standard array assignment 
        """
        #f.set_index(0, xdot) #dxdt
        f[0] = xdot
        if abs(V-xdot) > 0:
            f[1] = ((V-xdot)/abs(V-xdot)*0.5*self.C_D*rho*(V-xdot)**2*self.A_valve+deltap*self.A_valve-self.k_valve*x)/(self.m_eff) #d(xdot)dt
        else:
            f[1] = (deltap*self.A_valve-self.k_valve*x)/(self.m_eff) #d(xdot)dt
        return
        
    #@cython.cdivision(True)
    def _flux_dominant(self, f,x, xdot,rho,V):
        f[0] = xdot #dxdt
        if abs(V-xdot) > 0:
            f[1] = ((V-xdot)/abs(V-xdot)*0.5*self.C_D*rho*(V-xdot)**2*self.A_valve+(V-xdot)/abs(V-xdot)*rho*(V-xdot)**2*self.A_port-self.k_valve*x)/(self.m_eff) #d(xdot)dt
        else:
            f[1] = (-self.k_valve*x)/(self.m_eff) #d(xdot)dt
        return
    
    def set_xv(self, xv):
        self.xv = xv.copy()
        #If valve opening is less than zero, just use zero (the valve is closed)
        if self.xv[0] < -1e-15 and self.xv.get_index(1) < 1e-15:
            #print 'closed, desired position is',self.xv.get_index(0),' and velocity is',self.xv.get_index(1)
            self.xv[0] = 0.0
            self.xv[1] = 0.0
        #If it predicts a valve opening greater than max opening, just use the max opening
        elif self.xv[0] > self.x_stopper:
            self.xv[0] = self.x_stopper
            self.xv[1] = 0.0
        #TODO: stuck here
        print self.xv,'set_xv',xv
        
    def A(self):
        if self.xv is None:
            print 'self.xv is None'
        x = self.xv[0]
        if x >= self.x_tr:
            return self.A_port
        else:
            return pi*x*self.d_valve
    
    def flow_velocity(self, State_up, State_down):
        """
        For a given set of states, and a known valve lift, first
        check whether it is within the valve lift range, and then
        calculate the flow velocity
        """
        A = self.A()
        x = self.xv[0]
        if A > 0:
            if x > self.x_tr:
                return IsentropicNozzle(A, State_up, State_down, OUTPUT_VELOCITY)
            else:
                return x/self.x_tr*IsentropicNozzle(A, State_up, State_down, OUTPUT_VELOCITY)
        else:
            return 0.0
        
    #@cython.cdivision(True)
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
        #cdef double omega
        #TODO: check how to make it as normal python code
        #cdef arraym   
        f = empty_arraym(2)
        #cdef arraym
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
            
        print deltap, p_high, p_low, x, xdot, V, f,'valves'
            
        omega = Core.omega
        out_array[0] = f[0]/omega
        out_array[1] = f[1]/omega
        
        if abs(x) < 1e-15 and xdot < -1e-12:
            print 'stationary valve'
            out_array[0] = 0.0
            out_array[1] = 0.0
        
        return out_array #[dxdtheta, d(xdot)_dtheta]
    
    #TODO: check __cdict__
    def __dict__(self):
        items=['A_port','A_valve','d_valve','h_valve','d_port','m_eff','C_D','a_valve','l_valve',
               'rho_valve','k_valve','x_stopper','key_up','key_down','xv','State_up', 'State_down',
               'E']
        return {item:getattr(self,item) for item in items}
    
    def __repr__(self):
        """
        Return a representation of Valve Model for outputting to screen
        """
        s=''
        for item in self.__cdict__():
            s += item+' : '+str(getattr(self,item))+'\n'
        return s
            
    def __reduce__(self):
        return rebuildValveModel,(self.__cdict__(),) #TODO: check __cdict__
        
    
if __name__ == '__main__':
    
    print 'add function to be run'
    