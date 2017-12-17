from __future__ import division, print_function, absolute_import
from math import pi,cos,sin,sqrt
from PDSim.misc.scipylike import trapz
from PDSim.flow import flow_models
from PDSim.core.core import PDSimCore
from PDSim.misc.datatypes import arraym
from ._recip import _Recip

class Recip(PDSimCore,_Recip):
    """
    Recip is derived from :class:`PDSimCore <PDSim.core.core.PDSimCore>`.
    
    functions V_dV and heat_transfer_callback are provided in _Recip Cython class 
    """
    def __init__(self):
        PDSimCore.__init__(self)
        
    def __getstate__(self):
        """
        A function for preparing class instance for pickling
         
        Combine the dictionaries from the _Recip base class and the Recip
        class when pickling
        """
        py_dict = self.__dict__.copy()
        py_dict.update(_Recip.__cdict__(self))
        return py_dict

    def __setstate__(self, d):
        """
        A function for unpacking class instance for unpickling
        """
        for k,v in d.iteritems():
            setattr(self,k,v)
            
    def V_dV(self, theta):
        "A thin wrapper around Cython code for pickling purposes"
        return _Recip.V_dV(self, theta)
    
    def Vdisp(self):
        """
        Returns displacement of compressor in m\ :math:`^3`\ /revolution
        """
        return 2*self.crank_length*self.A_piston
    
    def V_shell(self, theta):
        
        return self.shell_volume, 0
        
    def TubeCode(self,Tube,**kwargs):
        Tube.Q = flow_models.IsothermalWallTube(Tube.mdot,Tube.State1,Tube.State2,Tube.fixed,Tube.L,Tube.ID,T_wall=self.Tlumps[0])
    
    def Inlet(self, FlowPath, **kwargs):
        try:
            FlowPath.A=0.01**2/4*pi
            mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
        
    def PistonLeakage(self, FlowPath, **kwargs):
        try:
            FlowPath.A=self.delta_gap*self.piston_diameter*pi
            mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
            return mdot
        except ZeroDivisionError:
            return 0.0
                
    def Suction(self,FlowPath,**kwargs):
        if FlowPath.key_up=='A':
            ## pressure in compressor higher than the inlet line
            ## valve is closed - no flow
            return 0.0
        else:
            try:
                FlowPath.A=self.suction_valve.A()
                #FlowPath.A=self.A_suction
                mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
                return mdot
            except ZeroDivisionError:
                return 0.0
        
    def Discharge(self,FlowPath,**kwargs):
        if FlowPath.key_down=='A':
            ## pressure in compressor lower than the discharge line
            ## valve is closed - no flow
            return 0.0
        else:
            try:
                FlowPath.A=self.discharge_valve.A()
#                FlowPath.A=self.A_discharge
                mdot=flow_models.IsentropicNozzle(FlowPath.A,FlowPath.State_up,FlowPath.State_down)
                
                return mdot
            except ZeroDivisionError:
                return 0.0

    def heat_transfer_callback(self,theta):
        "A thin wrapper of the Cython code for pickling purposes"
        return _Recip.heat_transfer_callback(self,theta)   
     
#    def heat_transfer_callback(self,theta):
#        T_w  = self.Tlumps[0] #[K]
#        V,dV=self.V_dV(theta) #[m3,m3/radian]
#        D_h = self.piston_diameter #[m]
#        A_ht = pi*self.piston_diameter*(V/self.A_piston) #[m2]
#        
#        Pr = self.CVs['A'].State.Prandtl #[-]
#        rho = self.CVs['A'].State.rho #[kg/m3]
#        k = self.CVs['A'].State.k #[kW/m-K]
#        mu = self.CVs['A'].State.visc #[Pa-s]
#        T = self.CVs['A'].State.T #[K]
#        u = abs(0.5*dV*self.omega/self.A_piston) #[m/s]
#        Re = (rho*u*D_h)/mu #[kg/m3*m/s*m/Pa/s]=[-]
#        
#        h_c = 0.053*(k/D_h)*Pr**(0.6)*Re**(0.8) #[kW/m2/K]
#        Q = h_c*A_ht*(T_w-T)   #Net heat into control volume [kW]
#        if self.CVs.N > 1:
#            return listm([Q,0])
#        else:
#            return listm([Q])
    
    def mechanical_losses(self):
        """
        
        Mathematical Analysis:
        
        .. math::
        
            x = {L_c}\cos \\theta  + \sqrt {L_1^2 - L_c^2{{\\left( {\\sin \\theta } \\right)}^2}}
            
        .. math::
        
            \\dot x =  - {L_c}\\sin \\theta \\dot \\theta  + \\frac{{ - L_c^2\\sin 2\\theta }}{{2\\sqrt {L_1^2 - L_c^2{{\\left( {\\sin \\theta } \\right)}^2}} }}\\dot \\theta
        
        .. math::
            
            \\bar \\dot x =  - \\frac{{\\dot \\theta }}{\\pi }\\int_0^\pi  {\\left[ { - {L_c}\\sin \\theta  + \\frac{{ - L_c^2\\sin 2\\theta }}{{2\\sqrt {L_1^2 - L_c^2{{\\left( {\\sin \\theta } \\right)}^2}} }}} \\right]d\\theta }
        
        .. math::
            
            \\bar \\dot x =  - \\frac{{\\dot \\theta }}{\\pi }\\left[ {{L_c}\\cos \\theta } \\right]_0^\pi  =  - \\frac{{{L_c}\\dot \\theta }}{\\pi }( - 1 - 1) = \\frac{{2{L_c}\\dot \\theta }}{\\pi }
        
        """
        #Oil with viscosity of 10 cSt (=10e-6 m^2/s) and density of 860 kg/m^3
        A_length = pi*(self.piston_diameter)*self.piston_length #[m2]
        u_ave = 2*self.crank_length/pi*self.omega
        F_viscous = (self.mu_oil*A_length*u_ave)/self.delta_gap #[N]
        W_dot_friction = (F_viscous*u_ave)/1000 #[kNm/s]=[kW]
        return W_dot_friction
        
    def ambient_heat_transfer(self):
        """
        The ambient heat transfer for the compressor in kW
        
        Returns a positive value if heat is added to the compressor
        """
        return self.h_shell*self.A_shell*(self.Tamb-self.Tlumps[0]) #[kW]
        
    def lump_energy_balance_callback(self):
        #Mean heat transfer rate from gas to lump using trapezoidal method for numerical integration
        #Note: Qdot_from_gas has opposite sign of heat transfer TO the gas as calculated in the heat_transfer_callback
        #Note: Qdot_from_gas will likely be less than 0 because heat is being removed and delivered to the gas
        self.Qdot_from_gas = -trapz(self.Q[0,0:self.Itheta],self.t[0:self.Itheta])/(self.t[self.Itheta-1]-self.t[0]) #[kW]
        #Mechanical losses are added to the lump
        self.Wdot_mechanical = self.mechanical_losses() #[kW]
        #Heat transfer between the shell and the ambient
        self.Qamb = self.ambient_heat_transfer() #[kW]
        #Note: heat transfer TO the gas in the tubes is given a positive sign
        #      so sign is flipped for the lump
        self.Qtubes = -sum([Tube.Q for Tube in self.Tubes])
        return self.Qdot_from_gas + self.Wdot_mechanical + self.Qamb + self.Qtubes
    
    def pre_solve(self):
        """
        Other calculations that are an indirect result of the inputs
        """
        self.A_piston = pi*(self.piston_diameter)**2/4
        self.V_dead=self.A_piston*self.x_TDC
        self.A_discharge=pi*self.d_discharge**2/4
        self.A_suction=pi*self.d_suction**2/4
        
    #Name gets mangled in the core base class, so un-mangle it
    def _PDSimCore__post_solve(self):
        """
        {\eta _{motor}} = \frac{{{{\dot W}_{shaft}}}}{{{{\dot W}_{shaft}} + {{\dot W}_{motor}}}}
        {\eta _{motor}}\left( {{{\dot W}_{shaft}} + {{\dot W}_{motor}}} \right) = {{\dot W}_{shaft}}
        {{\dot W}_{motor}} = \frac{{{{\dot W}_{shaft}}}}{{{\eta _{motor}}}} - {{\dot W}_{shaft}}
        """
        #Call the base class function
        PDSimCore._PDSimCore__post_solve(self)
        
        #Extra code for the recip
        #Motor losses
        self.Wdot_motor = self.Wdot*(1/self.eta_motor-1)
        
        #Electrical Power
        self.Wdot_electrical = self.Wdot + self.Wdot_motor
        
        #Overall isentropic efficiency
        self.eta_oi = self.Wdot_i/self.Wdot_electrical
        
if __name__=='__main__':
    print("Running this file doesn't do anything")
