from __future__ import division

cimport cython



"""
State_Flooded should receive   set_Fluid and set_Liquid as well as T,P from the 
refrigerant, meaning the State class of CoolProp

To add as inputs by the user:

__hasLiquid__ == True/False
set_Liquid()


"""

cdef class State_Flooded(State):
    
    def __init__(self, Ref = None,Liq = None,P,T,xL):
        
        self.Ref = Ref
        self.Liq = Liq
        self.P = P
        self.T = T
        self.xL = xL
    
    
    #Properties in PDSim - KSI & in CoolProp - SI
    ## Entropy        
    def s_gas(self,Ref,T,P):
        self.s_g = PropsSI('S','T',float(T),'P',float(P)*1000,Ref)  #J/Kg-K
        return self.s_g
    
    def s_liq(self,Liq,T):
        self.T0 = 273.15
        self.P0 = 101.325

        if Liq == "PAO":
            self.s_l = 1.940 * log(T/self.T0)
            
        elif Liq == "PAG":
            try:
                self.s_l = 2.74374E-03*(T-self.T0)+1.08646*log(T/self.T0)
            except:
                a=4
                
        elif Liq == "POE":
            self.s_l = 2.30 * log(T/self.T0)
            
        elif Liq == "Duratherm_LT":
            self.s_l = (3.4014*(T-298)+1094.3*log(T/298.0))/1000 
            
        elif Liq == "Water":
            cl_A=92.053
            cl_B=-0.039953
            cl_C=-0.00021103
            cl_D=5.3469E-07
            MM_l=18.0153
            try:
                self.s_l = (cl_A*log(T/298.15) + cl_B*(T-298.15) + cl_C/2.0*(T*T-298.15*298.15) + cl_D/3.0*(T*T*T-298.15*298.15*298.15))/MM_l
            except:
                a=1
                
        else:
            print "Invalid fluid"
    
        return self.s_l*1000  #J/kg-K
    
    def s_mix(self,Ref,Liq,T,P,xL):  #Entropy of the mixture as a function of temperature [K] and pressure [kPa].  Output in kJ/kg-K
        self.s_m = xL*self.s_liq(Liq,T) + (1-xL)*self.s_gas(Ref,T,P)
        
        if isnan(self.s_m) == True:
            print 's_m is a NaN'
        if isinf(self.s_m) == True:
            print 's_m is Infinite'
        else:
            return self.s_m
    
    ## Internal Energy
    def u_gas(self,Ref,T,P):
        
        self.u_g = PropsSI('U','T',float(T),'P',float(P)*1000,Ref)   #J/kgK
        return self.u_g
        
    def u_liq(self,Liq,T):
        
        if Liq == "Duratherm_LT":
            #internal energy [kJ/kg] of Duratherm LT given T in K
            self.u_l = (3.4014/2*pow(T,2)+1094.3*T)/1000  #LT
            
        elif Liq == "Zerol60":
            self.u_l = (5.186/2*pow(T,2)+337.116*T)/1000  #Zerol60
            
        elif Liq == "POE":
            self.u_l = (2.0935/2*T**2 + 1186.7*T)/1000 #POE 150 SUS
            
        else:
            print "Invalid Fluid"
        
        return self.u_l*1000  #J/kg
    
    def u_mix(self,Ref,Liq,T,P,xL):  #input in K, [-] output in J/kg
        
        self.u_l = self.u_liq(Liq,T)
        self.u_g = self.u_gas(Ref,T,P)
        
        self.u_m = xL*self.u_l + (1-xL)*self.u_g
        
        if isnan(self.u_m) == True:
            print 'u_m is a NaN'
        if isinf(self.u_m) == True:
            print 'u_m is Infinite'
        else:
            return self.u_m  #J/kg
    
    ## Enthalpy
    def h_gas(self,Ref,T,P):

        self.h_g = PropsSI('H','T',float64(T),'P',float64(P)*1000,Ref)   #J/kg
        return self.h_g
    
    def h_liq(self,Liq,T,P):  # Enthalpy of the mixture as a function of temperature [K] and pressure [kPa].  Output in kJ/kg
        self.T0 = 273.15
        self.P0 = 101.325
        self.h = 0
        self.h_l = 0
    
        if Liq =='PAO':
            self.h_l = 1.940*(T-self.T0)+(P-self.P0)/849
        elif Liq=='PAG':
            # PAG 0-OB-1020 from Tribology Data Handbook
            self.rho_l = -0.726923*float64(T)+1200.22;
            self.h_l = 2.74374E-03*(float64(T)**2-self.T0**2)/2.0+1.08646*(float64(T)-self.T0)+(float64(P)-self.P0)/self.rho_l;
        elif Liq == 'POE':
            # From Totten, p 261, cp=0.55 cal/g-C --> 2.30 kJ/kg-K
            self.h_l = 2.30*(T-self.T0)+(P-self.P0)/930
        elif Liq == 'Duratherm_LT':
            #the specific enthalpy of Duratherm LT [kJ/kg-k]"
            self.h_l = self.u_liq(Liq,T)/1000 + (P-self.P0)/self.rho_liq(Liq,T)
        elif Liq =='Water':
            cl_A = 92.053
            cl_B = -0.039953
            cl_C = -0.00021103
            cl_D = 5.3469E-07
            MM_l = 18.0153
            self.h_l = (cl_A*(T-298.15) + cl_B/2.0*(T**2-298.15**2) + cl_C/3.0*(T**3-298.15**3) + cl_D/4.0*(T**4-298.15**4))/MM_l+(P-self.P0)/rho_m(Ref,Liq,T,P,1.0)
        else:
            print "Invalid fluid"
            
        return self.h_l*1000
        
    def h_mix(self,Ref,Liq,T,P,xL):
        
        self.h_l = self.h_liq(Liq,T,P)
        self.h_g = self.h_gas(Ref,T,P)
        
        self.h_m = xL*self.h_l + (1-xL)*self.h_g
    
        if isnan(self.h_m) == True:
            print 'h_m is a NaN'
        if isinf(self.h_m) == True:
            print 'h_m is Infinite'
        else:
            return self.h_m
    
    ## Density
    def rho_gas(self,Ref,T,P):    
        self.rho_g = PropsSI('D','T',float(T),'P',float(P)*1000,Ref)    
        return self.rho_g
        
    def rho_liq(self,Liq,T):
        if Liq == 'PAO':
            self.rho_l = 849
        elif Liq == 'PAG':  # PAG 0-OB-1020 from Tribology Data Handbook
            self.rho_l = -0.726923*T+1200.22;
        elif Liq == 'POE':
            self.rho_l = 930
        elif Liq == 'Duratherm_LT':
        #density [kg/m^3] of Duratherm LT given T in K"
            self.rho_l = -0.6793*T + 1012.4 
        elif Liq =="Water":
            # Water Props from Yaws
            rhol_A=0.3471     
            rhol_B=0.274      
            rhol_n=0.28571    
            rhol_Tc=647.13
            self.rho_l = rhol_A/pow(rhol_B,pow(1-T/rhol_Tc,rhol_n))*1000;
        else:
            print "Invalid fluid"    
        return self.rho_l
    
    def rho_mix(self,Ref,Liq,T,P,xL,**kwargs):   #Density of the mixture assuming homogeneous mixture properties
    
        self.rho_l = self.rho_liq(Liq,T)
        self.rho_g = self.rho_gas(Ref,T,P)
        
        self.rho = None
        self.S = None
        self.alpha = None
        self.v_l = 1.0/self.rho_l
        self.v_g = 1.0/self.rho_g
        x_g = 1.0 - xL
        x = x_g
        
        if 'model' not in kwargs or kwargs['model']=='HEM':    
            self.S = 1
        elif kwargs['model']=='Zivi':
            self.S = (self.v_g/self.v_l)**(0.33333) #Eqn. 4.57 from Chisholm
        elif kwargs['model']=='Fauske':
            self.S = (self.v_g/self.v_l)**(0.5) #Eqn. 4.57 from Chisholm
        
        self.rho = (x + self.S*(1.0-x))/(x*self.v_g + self.S*(1.0-x)*self.v_l) #Eq 2.36 from Chisholm
        
        if x>0:
            self.alpha = 1.0/(1.0+(1.0-x)/x*self.rho_g/self.rho_l*self.S)
        else:
            self.alpha=0.0
        
        if 'alphaOn' in kwargs and kwargs['alphaOn']==True:
            return self.rho,self.alpha
        else:
            return self.rho
     
    ## Specific Heat        
    def cp_gas(self,Ref,T,P):
    
        self.cp_g = PropsSI('C','T',float(T),'P',float(P)*1000,Ref)   #J/kgK
        return self.cp_g

    def cp_liq(self,Liq,T):
    
        if Liq ==  'PAO':
            self.cp_l = 1.940
        elif Liq == 'POE':
            self.cp_l = 2.30
        elif Liq == 'PAG':
            # PAG 0-OB-1020 from Tribology Data Handbook
            # T in K, cp in kJ/kg-K
            self.cp_l = 2.74374E-03*T+1.08646;
        elif Liq == 'Duratherm_LT':
            #specific heat [kJ/kg-K] of Duratherm LT given T in K
            self.cp_l = (3.4014*T + 1094.3)/1000
        elif Liq == "Water":
            MM_l=18.0153
            cl_A=92.053/MM_l
            cl_B=-0.039953/MM_l
            cl_C=-0.00021103/MM_l
            cl_D=5.3469E-07/MM_l
            self.cp_l = cl_A + cl_B*T + cl_C*T*T + cl_D*T*T*T
        elif Liq == "ACD100FY":
            # 273 < T [K] < 387
            self.cp_l = 1.304 + 1.035e-3*T+2.801e-6*T**2   #[kJ/kg-K]
        return self.cp_l*1000
            
    def cp_mix(self,Ref,Liq,T,P,xL):
        
        xL = float(xL)
        self.cp_g = self.cp_gas(Ref,T,P)
        self.cp_l = self.cp_liq(Liq,T)
        
        self.cp_m = xL*self.cp_l + (1.0-xL)*self.cp_g
    
        if isnan(self.cp_m) == True:
            print 'cpm is a NaN'
        if isinf(self.cp_m) == True:
            print 'cpm is Infinite'
        else:
            return self.cp_m
            
    ##        
    def cv_gas(self,Ref,T,P):
        
        self.cv_g = PropsSI('O','T',float(T),'P',float(P)*1000,Ref)   #J/kgK
        return self.cv_g
    
    def cv_mix(self,Ref,Liq,T,P,xL):
        
        xL = float(xL)
        self.cv_g = self.cv_gas(Ref,T,P)
        self.cv_l = self.cp_liq(Liq,T)
        
        self.cv_m = xL*self.c_l + (1.0-xL)*self.cv_g
    
        if isnan(self.cv_m) == True:
            print 'cvm is a NaN'
        if isinf(self.cv_m) == True:
            print 'cvm is Infinite'
        else:
            return self.cv_m
    
    ## Dynamic Viscosity [Pa-s]
    def mu_liq(self,Liq,T):
        self.mu_l = 0
        
        if Liq == 'Duratherm_LT': 
            self.mu_l = 8e12*pow(T,-6.001)  #LT
    
        elif Liq == 'Zerol60': 
            self.mu_l = 1.0*(-0.0001235*T + 0.04808) #Zerol60 Vincent
        
        elif Liq == 'POE':
            #POE equation only valid from 60 C to 120 C
            #POE 150 SUS
            self.mu_l = 0.000000000517*T**4 - 0.000000795840*T**3 + 0.000460766590*T**2 - 0.118976538068*T + 11.571730524692                   
        elif Liq == "ACD100FY":
            # 313 < T [K] < 387
            #mu_l = 1.603e+46*pow(T,-19.69) #T [K]   
            self.mu_l = 2.0022e+31*pow(T,-1.2961e+01)
        else:
            print "Invalid fluid"
            
        return self.mu_l
        
    def mu_mix(self,Ref,Liq,T,P,xL):
        self.mu_l = self.mu_liq(Liq,T)
        self.mu_g = PropsSI('V','T',float(T),'P',float(P)*1000,Ref)
        self.mu_m = 1/(xL/self.mu_l + (1-xL)/self.mu_g )
        if isnan(self.mu_m) == True:
            print 'mu_m is a NaN'
        if isinf(self.mu_m) == True:
            print 'mu_m is Infinite'
        else:
            return 1/(xL/self.mu_l + (1-xL)/self.mu_g )
            
    ## Void Fraction
    def VoidFrac(self,Ref,Liq,T,P,xL):   #Input in K,kPa , [-] output [-]   
        self.rho_l = self.rho_liq(Liq,T)
        self.rho_g = self.rho_gas(Ref,T,P)
        self.VF = ((1-xL)/self.rho_g)/((1-xL)/self.rho_g + xL/self.rho_l)
        
        if isnan(self.VF) == True:
            print 'VF is a NaN'
        if isinf(self.VF) == True:
            print 'VF is Infinite'
        else:
            return self.VF
    
    ## Thermal Conductivity [W/m-K]
    def k_gas(self,Ref,T,P):
        self.k_g = PropsSI('L','T',float(T),'P',float(P)*1000,Ref)
        return self.k_g
    
    def k_liq(self,Liq,T):
    
        if Liq == 'Duratherm_LT': 
            self.k_l = -9e-5*T + 0.1223
    
        elif Liq == 'Zerol60': 
            self.k_l =  0.1153 #Zerol60 Vincent
        
        elif Liq == 'POE':  
            self.k_l = 0.138   #POE 150 SUS               
        
        elif Liq == "ACD100FY":
            self.k_l = 0.13    #TODO: check value of conductivity
        
        else:
            print "Invalid fluid"
            
        return self.k_l
        
    def k_mix(self,Ref,Liq,T,P,xL):  #Input in K, kPa ,[-] output in kW/m-K
        self.k_l = self.k_liq(Liq,T)
        self.k_g = self.k_gas(Ref,T,P)
        self.VF = self.VoidFrac(Ref,Liq,T,P,xL)
        self.k_m = (1-self.VF)*self.k_l + self.VF*self.k_g
        if isnan(self.k_m) == True:
            print 'k_m is a NaN'
        if isinf(self.k_m) == True:
            print 'k_m is Infinite'
        else:
            return (1-self.VF)*self.k_l + self.VF*self.k_g
    
    ## Prandtl Number
    def Pr_mix(self,Ref,Liq,T,P,xL): #TODO: maybe it is Prm (pag491 Ian Thesis)
        self.Pr_m = (self.cp_mix(Ref,Liq,T,P,xL)*self.mu_mix(Ref,Liq,T,P,xL))/self.k_mix(Ref,Liq,T,P,xL)
        
        if isnan(self.Pr_m) == True:
            print 'Pr_m is a NaN'
        if isinf(self.Pr_m) == True:
            print 'Pr_m is Infinite'
        else:
            return self.Pr_m
            
    ##
    def kstar_mix(self,Ref,Liq,T,P,xL):
    
        self.kstar_m =((1-xL)*self.cp_gas(Ref,T,P) + xL*self.cp_liq(Liq,T))/((1-xL)*self.cv_gas(Ref,T,P) + xL*self.cp_liq(Liq ,T))
        
        if isnan(self.kstar_m) == True:
            print 'kstar_m is a NaN'
        if isinf(self.kstar_m) == True:
            print 'kstar_m is Infinite'
        else:
            return self.kstar_m
            
    ##
    def e_mix(self,Ref,Liq,T,P,xL):
        
        self.T0 =25+273.15
        self.P0 =101.325
        
        self.e_m = (self.h_mix(Ref,Liq,T,P,xL) - self.h_mix(Ref,Liq,self.T0,self.P0,xL)) - self.T0*(self.s_mix(Ref,Liq,T,P,xL) - self.s_mix(Ref,Liq,self.T0,self.P0,xL));
    
        if isnan(self.e_m) == True:
            print 'e_m is a NaN'
        if isinf(self.e_m) == True:
            print 'e_m is Infinite'
        else:
            return float(self.e_m)    
        
    ##
    def dpdT_const_V(self,Ref,Liq,T,P1,xL):
        
        global rho_mix
        
        delta =1e-5
        self.v = 1/self.rho_mix(Ref,Liq,T,P1,xL)
        self.f = lambda P2: 1.0/self.rho_mix(Ref,Liq,T+delta,P2[0],xL) - self.v
        
        P2 = optimize.fsolve(self.f,P1)
        
        if isnan((P2-P1)/delta) == True:
            print 'dpdt_v is a NaN'
        if isinf((P2-P1)/delta) == True:
            print 'dpdt_v is Infinite'
        else:
            return float((P2-P1)/delta)   
    
    ##        
    def dudxL_mix(self,Ref,Liq,T,P,xL):

        delta =.001;
        self.dudxL_m = (self.u_mix(Ref,Liq,T,P,xL+delta) - self.u_mix(Ref,Liq,T,P,xL))/delta
        
        return self.dudxL_m
        
        
    def T_sp(self,Ref,Liq,s,p,xL,T_guess):   
    #Solve for the temperature which gives the same entropy - s [kJ/kg-K], p [kPa], T [K]
    
        global s_mix
        """
        When fsolve() is called, it will use a single element
        ndarray but we only want a double to be passed along,
        so take the first element of the 1-element array, which
        is a 64-bit float
        """ 
        self.f = lambda T: self.s_mix(Ref,Liq,T[0],p,xL) - s
        self.T = optimize.fsolve(self.f,T_guess)
        if isnan(self.T) == True:
            print 'T_sp is a NaN'
        if isinf(self.T) == True:
            print 'T_sp is Infinite'
        else:
            return float(self.T)
    
    def T_hp(self,Ref,Liq,h,p,xL,T_guess):   
    # Solve for the temperature which gives the same enthalpy - h [kJ/kg], p [kPa], T [K]
    
        global h_mix
        
        """
        When fsolve() is called, it will use a single element
        ndarray but we only want a double to be passed along,
        so take the first element of the 1-element array, which
        is a 64-bit float
        """
        
        self.f = lambda T: self.h_mix(Ref,Liq,T[0],p,xL) - h
        self.T = optimize.fsolve(self.f,T_guess)
        if isnan(self.T) == True:
            print 'T_hp is a NaN'
        if isinf(self.T) == True:
            print 'T_hp is Infinite'
        else:
            return float(self.T)
    
    
    def Q_Th(self,Ref,Liq,T,h,Q_guess):
        
        """
        Solve for the two-phase quality of refrigerant in the dome for a given
        enthalpy and saturation temperature
        Caveat: Function only good for pure refrigerant
        """
        """
        When fsolve() is called, it will use a single element
        ndarray but we only want a double to be passed along,
        so take the first element of the 1-element array, which
        is a 64-bit float
        """
        
        print 'h',h
        self.h_l = PropsSI('H','T',T,'Q',0,Ref)
        self.h_v = PropsSI('H','T',T,'Q',1,Ref)
        self.Qth = (h - self.h_l)/(self.h_v - self.h_l)
        return self.Qth
    
    def Q_Ts(self,Ref,Liq,T,s,Q_guess):
        
        """
        Solve for the two-phase quality of refrigerant in the dome for a given
        enthalpy and saturation temperature
        Caveat: Function only good for pure refrigerant
        """
        """
        When fsolve() is called, it will use a single element
        ndarray but we only want a double to be passed along,
        so take the first element of the 1-element array, which
        is a 64-bit float
        """
        
        self.f = lambda Q : PropsSI('S','T',T,'Q',Q[0],Ref) - s
        self.Q = optimize.fsolve(self.f,Q_guess)
        return float(self.Q)
    
    
    def T_crit(self,Ref):    # Critical Temperature of the refrigerant [K] - From coolprop  ->  Props(Fluid,PropName)
        
        if Ref=='R245fa':
            self.T = PropsSI('R245fa', 'Tcrit')
        elif Ref=='R410A':
            self.T = PropsSI('R410A', 'Tcrit')
        elif Ref=='R744':
            self.T = 304.128
        elif Ref=='R404A':
            self.T = PropsSI('R404A', 'Tcrit')
        elif Ref=='R134a':
            self.T = PropsSI('R134a', 'Tcrit')      #Tcr=374.2
        elif Ref=='SES36':
            self.T = PropsSI('SES36', 'Tcrit')
        elif Ref=='R290':
            self.T = 369.82
        elif Ref=='R717':
            self.T = 405.4
        else:
            print 'Uh oh...  Refrigerant not found T_crit'
        return self.T
    
    def p_crit(self,Ref):    # Critical pressure of the refrigerant [kPa] - From coolprop  ->  Props(Fluid,PropName)
        if Ref=='R245fa':
            self.p = PropsSI('R245fa', 'pcrit')
        elif Ref=='R410A':
            self.p = PropsSI('R410A', 'pcrit')
        elif Ref=='R744':
            self.p = 304.128
        elif Ref=='R404A':
            self.p = PropsSI('R404A', 'pcrit')
        elif Ref=='R134a':
            self.p = PropsSI('R134a', 'pcrit')      #Tcr=374.2
        elif Ref=='SES36':
            self.p = PropsSI('SES36', 'pcrit')
        elif Ref=='R290':
            self.p = 369.82
        elif Ref=='R717':
            self.p = 405.4
        else:
            print 'Uh oh...  Refrigerant not found T_crit'
        return self.p
        
        
"""
#In Cython we would do something equivalent to CoolProp.pxd 
# This class inheritate from State_Flooded which calculates the mixture properties 
# and returns the values in a similar way of State class in CoolProp.pxd

cdef class StateFlood:
    cdef AbstractState pAS
    cdef readonly bytes Fluid, phase
    cdef int iFluid,iParam1,iParam2,iOutput
    cdef int iLiquid    # Maybe ??
    
    cdef double T_, rho_, p_, xL_,
    
    cpdef set_Fluid(self, string Fluid, string backend)
    cpdef set_Liquid(self, string Fluid)  #No need for backend
    
    #Now the properties are all defined above and we can be consistend with the definition
    #I used in containers_20150419  @properties
    
    cpdef speed_test(self, int N)
    cpdef update(self, dict params)
    cpdef update_ph(self, double p, double h)
    cpdef update_Trho(self, double T, double rho)
    cpdef State copy(self)
    cpdef double Props(self, constants_header.parameters iOutput) except *
    cpdef long Phase(self) except *
    cpdef double get_Q(self) except *
    cpdef double get_T(self) except *
    cpdef double get_p(self) except *
    cpdef double get_h_m(self) except *
    cpdef double get_rho_m(self) except *
    cpdef double get_s_m(self) except *
    cpdef double get_u_m(self) except *
    cpdef double get_visc_m(self) except *
    cpdef double get_cond_m(self) except *
    cpdef double get_cp_m(self) except *
    cpdef double get_cp0(self) except *   # TO CHECK
    cpdef double get_cv_m(self) except *
    cpdef double get_dpdT(self) except *



"""


