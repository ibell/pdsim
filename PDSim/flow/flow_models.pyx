
#cython: embedsignature=True

from __future__ import division
from CoolProp.State import State

import numpy as np
import cython
cimport cython


#Uncomment this line to use the python math functions
#from math import log,pi,e,pow,sqrt

"""
This is the module that contains all the flow models
"""

cdef public enum:
    TYPE_RADIAL
    TYPE_FLANK

cdef double ar_0 = 25932.1070099
cdef double ar_1 = 0.914825434095
cdef double ar_2 = -177.588568125
cdef double ar_3 = -0.237052788124
cdef double ar_4 = -172347.610527
cdef double ar_5 = -12.0687599808
cdef double ar_6 = -0.0128861161041
cdef double ar_7 = -151.202604262
cdef double ar_8 = -0.999674457769
cdef double ar_9 = 0.0161435039267
cdef double ar_10 = 0.825533456725
cdef double Re_star_radial = 5243.58194594

cdef double af_0 = -2.63970395918
cdef double af_1 = -0.567164431229
cdef double af_2 = 0.83655499929
cdef double af_3 = 0.810567167521
cdef double af_4 = 6174.02825667
cdef double af_5 = -7.60907962464
cdef double af_6 = -0.510200923053
cdef double af_7 = -1205.17482697
cdef double af_8 = -1.02938914174
cdef double af_9 = 0.689497785772
cdef double af_10 = 1.09607735134
cdef double Re_star_flank = 826.167177885
    
cdef class PyFlowFunctionWrapper(FlowFunctionWrapper):
    """
    This function is defined at the in order to do high-level 
    mass flow functions that are slower but can still be accessed by the 
    Cython code. 
    """
    def __init__(self,Function,kwargs):
        self.Function = Function
        self.kwargs = kwargs
        
    cpdef double call(self,FlowPath FP):
        return self.Function(FP,**self.kwargs)
    
    def __reduce__(self):
        if not isinstance(self.Function,str):
            Function_str = self.Function.__name__
        else:
            Function_str = self.Function
        return makePyFlowFunctionWrapper, ({'kwargs': self.kwargs,'Function':Function_str},)

def makePyFlowFunctionWrapper(kwds):
    return PyFlowFunctionWrapper(**kwds)
        
cdef class IsentropicNozzleWrapper(FlowFunctionWrapper):
    """
    A wrapper that can be added to call the isentropic nozzle model
    if the flow area is constant
    """
    
    cpdef double call(self, FlowPath FP):
        """
        Returns the mass flow rate from the isentropic nozzle model
        """
        return IsentropicNozzle(FP.A,
                                FP.State_up,
                                FP.State_down)
    

cdef class FlowFunctionWrapper(object):
    """
    A wrapper to contain the function that will be called
    
    Two methods are provided, call(FP) and __call__(FP).  The special method
     __call__(FP) allows an instance of FlowFunctionWrapper to be called 
     directly like:: 
     
         FFW = FlowFunctionWrapper()
         FFW(FP)
         
     or you can call the call() method like::
     
         FFW.call(FP)
         
     Having both methods allows cython functions to stay at the C++ layer and not
     need to come back to python layer for speed 
     
    """
    cpdef double call(self, FlowPath FP):
        pass
        
    def __call__(self, FlowPath FP):
        """
        This special method calls the call() function of the derived class
        """
        return self.call(FP)
    
cpdef double pow(double x, double y):
    return x**y

cpdef IsothermalWallTube(mdot,State1,State2,fixed,L,ID,OD=None,HTModel='Twall',Tmean=None,T_wall=None,Q_add = 0.0,alpha=None):
    """
    In this tube model, one of the nodes is known (fixed), but the other is calculated based on heat transfer and pressure drop for a given mass flow rate
    
    Assumes flow is turbulent and fully developed.  The flow is almost certainly fully turbulent, but it is unlikely to be truly fully developed
    
    Parameters
    ----------
    mdot : float
        mass flow rate [kg/s]
    State1 : :class:`State <CoolProp.State.State>` instance
        State number 1 for the Tube
    State2 : :class:`State <CoolProp.State.State>` instance
        State number 2 for the Tube
    fixed : int, ``1`` or ``2``
        Which node is fixed
    L : float
        Length of tube [m]
    ID : float 
        Internal diameter of tube [m]
    OD : float, optional
        Outer diameter of tube [m], not used if fixed wall temperature
    HTModel : `'Twall'` or [other models] , optional
        Key for what type of model is used for the heat transfer
    Tmean : float, optional
        Mean fluid temperature for evaluation of fluid properties [K]
    T_wall : float, optional
        Temperature of wall [K]
    Q_add : float, optional
        Additional amount of heat that will be added to the fluid in the tube [kW]
        This term is not added to the amount of heat transfer returned from this function
    alpha : float, optional
        The heat transfer coefficient [kW/m2/K].  If not provided, calculated from correlations
        
    Returns
    -------
    Q : float
        The amount of heat transfer [W], not including ``Q_add``
        
    """
        
    #Use the provided value for Tmean if it is a float or integer
    if not (isinstance(Tmean,float) or isinstance(Tmean,int)):
        if fixed==1:
            Fluid=State1.Fluid
            Tmean=State1.T
        elif fixed==2:
            Fluid=State2.Fluid
            Tmean=State2.T
        else:
            print 'fixed not provided'
            #raise AttributeError
        
    if State1.hasLiquid:
        raise ValueError('Fill in here')
#        Pr=Pr_mix(Ref,Liq,T_in,p_in,xL_in) #[-]
#        Re=4.0*mdot/(PI*mu_mix(Ref,Liq,T_in,p_in,xL_in)*D) #[-]
#        k=k_mix(Ref,Liq,T_in,p_in,xL_in) #[kW/m-K]
#        cp=cp_mix(Ref,Liq,T_in,p_in,xL_in) #[kJ/kg-K]
#        hc=0.023*k/D*pow(Re,0.8)*pow(Pr,0.4) #[kW/m^2-K]
    else:
        if fixed==1:
            p=State1.p
            T1=State1.T
        elif fixed==2:
            p=State2.p
            T2=State2.T
            
        #Q_add needs to be in W
        Q_add *= 1000
        
        
        
        S=State(Fluid,{'T':Tmean,'P':p})
            
        mu = S.visc  #kg/m-s
        cp = S.cp*1000. #J/kg-K
        k = S.k*1000. #W/m-K
        rho = S.rho #kg/m^3
        Pr = cp * mu / k #[-]

        InnerFlowArea=pi*(ID**2)/4.0
        u=mdot/(InnerFlowArea*rho)
        Re=rho*u*ID/mu

        # Friction factor of Churchill (Darcy Friction factor where f_laminar=64/Re)
        e_D = 0
        A = ((-2.457 * log( (7.0 / Re)**(0.9) + 0.27 * e_D)))**16
        B = (37530.0 / Re)**16
        f = 8 * ((8/Re)**12.0 + 1 / (A + B)**(1.5))**(1/12)

        # Heat Transfer coefficient of Gnielinski
        Nu = (f/8)*(Re-1000)*Pr/(1+12.7*sqrt(f/8)*(Pr**(0.66666)-1)) #[-]
        
        # alpha is not provided directly, use the correlation
        if alpha is None:
            alpha = k*Nu/ID #W/m^2-K

        #Pressure gradient using Darcy friction factor
        G=mdot/InnerFlowArea
        dp_dz=-f/rho*G**2/(2*ID)
        DELTAP=dp_dz*L

        if fixed==1:
            #The outlet temperature considering just the wall heat transfer 
            T2_star=T_wall-(T_wall-T1)*exp(-pi*ID*L*alpha/(mdot*cp))
            
            #Get the actual outlet temperature based on the additional heat input
            T2 = T2_star + Q_add/(mdot*cp)
            
            # Q is defined to be positive if heat transferred from wall to fluid
            #
            # It only includes the term from the wall heat transfer 
            Q=mdot*cp*(T2_star-T1)
            
            State2.update({'T':T2,'P':p+DELTAP/1000})
        else:
            #Get the wall heat transfer outlet temperature based on the additional heat input
            T2_star = T2 - Q_add/(mdot*cp)
            
            T1=T_wall-(T_wall-T2_star)/exp(-pi*ID*L*alpha/(mdot*cp))
            
            # Q is defined to be positive if heat transferred from wall to fluid 
            Q=mdot*cp*(T2_star-T1)
            
            State1.update({'T':T1,'P':p-DELTAP/1000})
        
#        print 'T1',T1
#        print 'T2_star',T2_star
#        print 'T2',T2
#        print 'p',p
#        print 'mdot',mdot
#        print 'others: ID,L,alpha,cp',ID,L,alpha,cp
#        print 'DELTAP',DELTAP
#        print 'Fluid',Fluid
#        print 'Q_add', Q_add
#        print 'Q', Q
        
        return Q/1000.0

def rebuildValveModel(d):
    VM = ValveModel(d.pop('d_valve'),d.pop('d_port'),d.pop('C_D'),d.pop('h_valve'),
                          d.pop('a_valve'),d.pop('l_valve'),d.pop('rho_valve'),
                          d.pop('E'),d.pop('x_stopper'),
                          d.pop('key_up'),d.pop('key_down'))
    for item in d:
        setattr(VM,item,d[item])
    return VM
        
cdef class ValveModel(object):
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
    
    #Give Cython type definitions for terms
    @cython.locals(d_valve=cython.double,d_port=cython.double,C_D=cython.double,
                   h_valve=cython.double,a_valve=cython.double, l_valve=cython.double, 
                   rho_valve=cython.double,E=cython.double,x_stopper=cython.double,
                   key_up=cython.bytes,key_down=cython.bytes)
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
        
    cpdef get_States(self, Core):
        """
        Core is the main model core, it contains information that 
        is needed for the flow models
        """
        exists_keys=Core.CVs.exists_keys
        Tubes_Nodes=Core.Tubes.Nodes
        for key,Statevar in [(self.key_up,'State_up'),(self.key_down,'State_down')]:
            ## Update the pointers to the states for the ends of the flow path
            if key in exists_keys:
                setattr(self,Statevar,Core.CVs[key].State)
            elif key in Tubes_Nodes:
                setattr(self,Statevar,Tubes_Nodes[key])                  
    
    cpdef tuple _pressure_dominant(self,double x, double xdot, double rho, double V, double deltap):
        dxdt = xdot
        dxdot_dt = (0.5*self.C_D*rho*(V-xdot)**2*self.A_valve+deltap*self.A_valve-self.k_valve*x)/(self.m_eff)
        return (dxdt,dxdot_dt)
        
    cpdef tuple _flux_dominant(self,double x, double xdot, double rho, double V):
        dxdt = xdot
        dxdot_dt = (0.5*self.C_D*rho*(V-xdot)**2*self.A_valve+rho*(V-xdot)**2*self.A_port-self.k_valve*x)/(self.m_eff)
        return (dxdt,dxdot_dt)
    
    cpdef set_xv(self, list xv):
        #If xv[0] is less than zero, just use zero
        if xv[0]<0.0:
            xv[0]=0.0
            xv[1]=0.0
        #If it predicts a valve opening greater than max opening, just use the max opening
        elif xv[0]>self.x_stopper:
            xv[0]=self.x_stopper
            xv[1]=0.0
        self.xv=xv
        
    cpdef double A(self):
        return pi*self.xv[0]*self.d_valve
    
    cpdef double flow_velocity(self,State State_up, State State_down):
        """
        For a given set of states, and a known valve lift, first
        check whether it is within the valve lift range, and then
        """
        if State_up.get_p()<State_down.get_p():
            return 0.0
        try:
            #Need a dummy value for the area to get a flow velocity even when the valve is closed
            A_dummy=0.001 
            mdot,others=IsentropicNozzle(A_dummy,State_up,State_down,full_output=True)
            return others['v']
        except ZeroDivisionError:
            return 0.0
        
    @cython.cdivision(True)
    cpdef list derivs(self,Core): 
        x=self.xv[0]
        xdot=self.xv[1]
        self.get_States(Core)
        rho=self.State_up.get_rho()
        p_high=self.State_up.get_p()
        p_low=self.State_down.get_p()
        deltap=(p_high-p_low)*1000
        if deltap<0.0:
            return [0.0,0.0]
        V = self.flow_velocity(self.State_up, self.State_down)
        x_tr = 0.25*(self.d_port**2/self.d_valve)
        if x>=x_tr:
            f=self._pressure_dominant(x,xdot,rho,V,deltap)
        else:
            f=self._flux_dominant(x,xdot,rho,V)
        return [_dummy/Core.omega for _dummy in f] #[dxdtheta, dxdot_dtheta]
    
    cpdef dict __cdict__(self):
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
        return rebuildValveModel,(self.__cdict__(),)

@cython.cdivision(True)        
cpdef double IsentropicNozzle(double A, State State_up, State State_down, bytes other_output = None):
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
    other_output : string
        If 'v', returns the velocity
        
    Returns
    -------
    
    """
    
    #some cython type declarations
    cython.declare(T_up = cython.double, T_down = cython.double, mdot = cython.double)
    # Since ideal, R=cp-cv, and k=cp/cv
    cp = State_up.get_cp0()
    R = 8.314472/State_up.get_MM()*1000 #[J/kg/K]
    cv = cp-R/1000.0 #[kJ/kg/K]
    k = cp / cv
    
    
    p_up=State_up.get_p()
    T_up=State_up.get_T()
    p_down=State_down.get_p()
    
    # Speed of sound
    c=(k*R*T_up)**0.5
    # Upstream density
    rho_up=p_up*1000.0/(R*T_up)
    pr=p_down/p_up
    pr_crit=(1+(k-1)/2)**(k/(1-k))
    
    if pr > pr_crit: 
        # Mass flow rate if not choked [kg/s]
        mdot=A*p_up*1000.0/(R*T_up)**0.5*(2*k/(k-1.0)*pr**(2.0/k)*(1-pr**((k-1.0)/k)))**0.5
        # Throat temperature [K]
        T_down=T_up*(p_down/p_up)**((k-1.0)/k)
        # Throat density [kg/m3]
        rho_down=p_down*1000.0/(R*T_down)
        # Velocity at throat
        v=mdot/(rho_down*A)
        # Mach number
        Ma=v/c
    else:
        # Mass flow rate if choked
        mdot=A*rho_up*(k*R*T_up)**0.5*(1.+(k-1.)/2.)**((1+k)/(2*(1-k)))
        # Velocity at throat
        v=c
        # Mach Number
        Ma=1.0
    

    if other_output is None:
        return mdot
    elif other_output == 'v':
        return v
    elif other_output == 'Ma':
        return Ma
     
@cython.cdivision(True)
cpdef double FrictionCorrectedIsentropicNozzle(double A, State State_up, State State_down, double delta, int Type, double t = -1.0, double ro = -1.0, bint full_output = False):
    """
    Frictionally-corrected nozzle model - the so-called hybrid leakage model
    
    Parameters
    ----------
    A : float
        Flow area at the minimum area [\ :math:`m^2`\ ], equal to :math:`\delta_{flank}h` for flank, 
        and :math:`s_r\delta_{radial}` for radial leakage.  Intended for scroll
        compressor, but can be used with other compressors
    State_up : :class:`State <CoolProp.State.State>` instance
        The State instance corresponding to the upstream side of the flow path
    State_down : :class:`State <CoolProp.State.State>` instance
        The State instance corresponding to the downstream side of the flow path
    delta : float
        Gap width in meters
    Type : int
        One of flow_models.TYPE_RADIAL or flow_models.TYPE_FLANK
    t : float
        Scroll wrap thickness in m
    ro : float
        Orbiting radius in m
    full_output : boolean
        If False, just return float value of mdot. If True, return tuple of 
        mdot and others, where others is a dictionary of other outputs with the
        keys `` 
        
    Notes
    -----
    If Type is 'radial', t must be provided
    
    If Type is 'flank', ro must be provided
        
    Implements the frictionally-corrected method of 
    Bell et al., *A Computationally Efficient Hybrid Leakage Model for Modeling Leakage in Positive Displacement Compressors*, "2012 International Compressor Engineering Conference at Purdue University
    
    Journal publication to follow
    """
    
    cdef double mdot,Re,v,mdot_ratio
    #Get the flow velocity and mass flow rate using the Isentropic nozzle model
    mdot = IsentropicNozzle(A, State_up, State_down)
    
    if abs(mdot) < 1e-12:
        return mdot
    # Hydraulic diameter
    Dh=2 * delta
    # Viscosity for Re
    mu = State_up.get_visc()
    rho_up = State_up.get_rho()
    # Reynolds number
    v = mdot/rho_up/A
    Re=rho_up*v*Dh/mu

    if (Type == <int>TYPE_RADIAL and t<=0):
        raise ValueError("Type 'radial' provided, but thickness of scroll [{0:g}] is not positive".format(t))
    elif (Type == <int>TYPE_RADIAL and Re > 1e-12):
        Re_star=Re_star_radial
        xi=1.0/(1.0+exp(-0.01*(Re-Re_star)))
        Lstar=t/0.005
        delta_star=delta/10e-6
        mdot_ratio=ar_0*pow(Lstar,ar_1)/(ar_2*delta_star+ar_3)*(xi*(ar_4*pow(Re,ar_5)+ar_6)+(1-xi)*(ar_7*pow(Re,ar_8)+ar_9))+ar_10
    elif (Type == <int>TYPE_FLANK and ro <= 0):
        raise ValueError("Type 'flank' provided, but orbiting radius of scroll [{0:g}] is not positive".format(ro))
    elif (Type == <int>TYPE_FLANK and Re > 1e-12):
        Re_star=Re_star_flank
        xi=1.0/(1.0+exp(-0.01*(Re-Re_star)))
        Lstar=ro/0.005
        delta_star=delta/10e-6
        mdot_ratio=af_0*pow(Lstar,af_1)/(af_2*delta_star+af_3)*(xi*(af_4*pow(Re,af_5)+af_6)+(1-xi)*(af_7*pow(Re,af_8)+af_9))+af_10
        #print Re,xi,mdot_ratio,af_4*pow(Re,af_5)+af_6,af_7*pow(Re,af_8)+af_9
    else:
        mdot_ratio=1.0
    mdot=mdot/mdot_ratio
    
    
    return mdot
