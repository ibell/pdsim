from __future__ import division
from CoolProp.State import State
from math import log,exp,pi
import numpy as np
import cython
"""
This is the module that contains all the flow models
"""

def pow(a,b):
    # This function is provided so that Cython can properly
    # optimize calls for the powers
    return a**b
    
def sqrt(x):
    # This function is provided so that Cython can properly
    # optimize calls for the sqrt function
    return x**0.5
        
def IsothermalWallTube(mdot,State1,State2,fixed,L,ID,OD=None,HTModel='Twall',Tmean=None,T_wall=None):
    """
    In this tube model, one of the nodes is known (fixed), but the other is calculated based on heat transfer and pressure drop for a given mass flow rate
    
    Assumes flow is turbulent and fully developed.  The flow is almost certainly fully turbulent, but it is unlikely to be truly fully developed
    
    Parameters
    ----------
    mdot : float
        mass flow rate [kg/s]
    State1 : State instance
        State number 1 for the Tube
    State2 : State instance
        State number 2 for the Tube
    fixed : int, `1` or `2`
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
            
        print Tmean,p,Fluid
        
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
        alpha = k*Nu/ID #W/m^2-K

        #Pressure gradient using Darcy friction factor
        G=mdot/InnerFlowArea
        dp_dz=-f/rho*G**2/(2*ID)
        DELTAP=dp_dz*L

        if fixed==1:
            T2=T_wall-(T_wall-T1)*exp(-pi*ID*L*alpha/(mdot*cp))
            
            # Q is defined to be positive if heat transferred from wall to fluid 
            Q=mdot*cp*(T2-T1)
            
            State2.update({'T':T2,'P':p+DELTAP/1000})
        else:
            
            T1=T_wall-(T_wall-T2)/exp(-pi*ID*L*alpha/(mdot*cp))
            
            # Q is defined to be positive if heat transferred from wall to fluid 
            Q=mdot*cp*(T2-T1)
            
            State1.update({'T':T1,'P':p-DELTAP/1000})

        return Q/1000.0

def rebuildValveModel(d):
    VM = ValveModel(d.pop('d_valve'),d.pop('d_port'),d.pop('C_D'),d.pop('h_valve'),
                          d.pop('a_valve'),d.pop('l_valve'),d.pop('rho_valve'),
                          d.pop('E'),d.pop('x_stopper'),
                          d.pop('key_up'),d.pop('key_down'))
    for item in d:
        setattr(VM,item,d[item])
    return VM
        
class ValveModel(object):
    
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
        
    def get_States(self, Core):
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
    
    def _pressure_dominant(self,x,xdot,rho,V,deltap):
        dxdt = xdot
        dxdot_dt = (0.5*self.C_D*rho*(V-xdot)**2*self.A_valve+deltap*self.A_valve-self.k_valve*x)/(self.m_eff)
        return (dxdt,dxdot_dt)
        
    def _flux_dominant(self,x,xdot,rho,V):
        dxdt = xdot
        dxdot_dt = (0.5*self.C_D*rho*(V-xdot)**2*self.A_valve+rho*(V-xdot)**2*self.A_port-self.k_valve*x)/(self.m_eff)
        return (dxdt,dxdot_dt)
    
    def set_xv(self,xv):
        #If xv[0] is less than zero, just use zero
        if xv[0]<0.0:
            xv[0]=0.0
            xv[1]=0.0
        #If it predicts a valve opening greater than max opening, just use the max opening
        elif xv[0]>self.x_stopper:
            xv[0]=self.x_stopper
            xv[1]=0.0
        self.xv=xv
        
    def A(self):
        return pi*self.xv[0]*self.d_valve
    
    def flow_velocity(self, State_up, State_down):
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
    def derivs(self,Core): 
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
    
    def __cdict__(self):
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
def IsentropicNozzle(A,State_up,State_down,full_output=False):   
    """
    The mass flow rate is calculated by using isentropic flow model
    
    Parameters
    ----------
    
    A : double
        Throat area of the nozzle [m^2]
    State_up : ``State`` instance
        Upstream ``State`` instance
    State_down : ``State`` instance
        Downstream ``State`` instance
    full_output : bool
        Outputs just mdot (full_output=False), or mdot,others if full_output=True
        
    Returns
    -------
    
    """
    # Since ideal, R=cp-cv, and k=cp/cv
    cp=State_up.get_cp0()
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
        T_down=T_up*pow(p_down/p_up,(k-1.)/k)
        # Throat density [kg/m3]
        rho_down=p_down*1000.0/(R*T_down)
        # Velocity at throat
        v=mdot/(rho_down*A)
        # Mach number
        Ma=v/c
    else:
        # Mass flow rate if choked
        mdot=A*rho_up*(k*R*T_up)**0.5*pow(1.+(k-1.)/2.,(1+k)/(2*(1-k)))
        # Velocity at throat
        v=c
        # Mach Number
        Ma=1.0
    
    if full_output==True:
        otherparameters={'Ma':Ma,'v':v}
        return mdot,otherparameters
    else:
        return mdot
    
#    // Hydraulic diameter
#        Dh=2.*delta;
#    // Viscosity for Re
#        mu=mu_g(Gas,T_up,p_up);
#    // Reynolds number
#        *Re=rho_up*v*Dh/mu;
#
#    if (Type==CORRECTED_RADIAL_NOZZLE && *Re>1e-12)
#    {
#        a=a_radial;
#        Re_star=Re_star_radial;
#        xi=1.0/(1.0+exp(-0.01*(*Re-Re_star)));
#        Lstar=t/0.005;
#        delta_star=delta/10e-6;
#        mdot_ratio=a[0]*pow(Lstar,a[1])/(a[2]*delta_star+a[3])*(xi*(a[4]*pow(*Re,a[5])+a[6])+(1-xi)*(a[7]*pow(*Re,a[8])+a[9]))+a[10];
#    }
#    else if (Type==CORRECTED_FLANK_NOZZLE  && *Re>1e-12)
#    {
#        a=a_flank;
#        Re_star=Re_star_flank;
#        xi=1.0/(1.0+exp(-0.01*(*Re-Re_star)));
#        Lstar=ro/0.005;
#        delta_star=delta/10e-6;
#        mdot_ratio=a[0]*pow(Lstar,a[1])/(a[2]*delta_star+a[3])*(xi*(a[4]*pow(*Re,a[5])+a[6])+(1-xi)*(a[7]*pow(*Re,a[8])+a[9]))+a[10];
#    }
#    else
#    {
#        mdot_ratio=1.0;
#    }
#    *mdot=*mdot/mdot_ratio;
#}
