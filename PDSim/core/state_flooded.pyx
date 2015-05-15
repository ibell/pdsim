from __future__ import division

cimport cython
import CoolProp.CoolProp
import math
from math import floor
from math import log
import warnings
from scipy import optimize
from numpy import float64,isnan,isinf

from CoolProp.CoolProp cimport constants_header

cdef dict paras = {constants_header.iDmass : 'D',
                   constants_header.iQ : 'Q',
                   constants_header.imolar_mass : 'M',
                   constants_header.iT : 'T',
                   constants_header.iHmass : 'H',
                   constants_header.iP : 'P',
                   constants_header.iCpmass : 'C',
                   constants_header.iCp0mass : 'C0',
                   constants_header.iCvmass : 'O',
                   constants_header.iviscosity : 'V',
                   constants_header.iconductivity : 'L',
                   constants_header.ispeed_sound: 'A',
                   constants_header.iSmass : 'S',
                   constants_header.iUmass : 'U'
}

cdef dict paras_inverse = {v:k for k,v in paras.iteritems()}

cdef class StateFlooded(State):
    
    def __init__(self, object Ref, string Liq, double P, double T, double xL, string model):
        cdef string backend = b'?'
        self.Fluid = Ref # Defined in base class
        self.p_ = P
        self.T_ = T
        
        self.Liq = Liq
        self.xL_ = xL
    
        self.model = model
        
        if b'::' in Ref:
            backend, Ref = Ref.split(b'::', 1)

        self.set_Fluid(Ref, backend)

    #cpdef update_with_guesses(self, constants_header.input_pairs ipair, double Value1, double Value2, double Value3, PyGuessesStructure guesses):
    #    """ Update function - wrapper of c++ function :cpapi:`CoolProp::AbstractState::update` """
    #    cdef cAbstractState.GuessesStructure _guesses
    #    _guesses.T = guesses.T
    #    _guesses.p = guesses.p
    #    _guesses.rhomolar_liq = guesses.rhomolar_liq
    #    _guesses.rhomolar_vap = guesses.rhomolar_vap
    #    _guesses.x = guesses.x
    #    _guesses.y = guesses.y
    #    self.thisptr.update_with_guesses(ipair, Value1, Value2, Value3, _guesses)

    ##
    cpdef update_TrhoxL(self, double T, double rho, double xL):
        """
        Just use the temperature, density and liquid fraction directly for speed

        Parameters
        ----------
        T: float
            Temperature [K]
        rho: float
            Mixture Density [kg/m^3]
        xL: float
            Liquid mass fraction [-]
        """
        self.T_ = T
        self.rho_ = rho
        self.xL_ = xL
        #self.thisptr.update(DmassTxL_INPUTS, rho, T, xL)
    
    ##
    cpdef update(self, dict params):
        """
        Parameters
        params, dictionary
            A dictionary of terms to be updated, with keys equal to single-char inputs to the Props function,
            for instance ``dict(T=298, P = 101.325)`` would be one standard atmosphere
        """
       
        # Convert to integer_pair input

        cdef double p, val1, val2, o1 = 0, o2 = 0
        cdef long iInput1, iInput2
        cdef bytes errstr
        cdef constants_header.input_pairs input_pair
        
        # What you need to do in this function is
        # a) For a given set of inputs, convert the triple of inputs into T (same for both phases), rho_g, rho_l, xL
        # b) Update the pAS using something like:
        #      self.pAS.update(constants_header.DmassT_INPUTS, rho_g, T)
        # c) Figure out the mixture T, p
        # d) Cache the mixture temperature and pressure in the variables self.p_ and self.T_ (defined in base class)

        # Convert inputs to input pair
        #items = list(params.items())
        #key1 = paras_inverse[items[0][0]]
        #key2 = paras_inverse[items[1][0]]
        
        # Convert to SI units
        #val1 = toSI(key1, items[0][1])
        #val2 = toSI(key2, items[1][1])
        
        #input_pair = _generate_update_pair(key1, val1, key2, val2, o1, o2)
        #self.pAS.update(input_pair, o1, o2);

        #self.T_ = self.pAS.T()
        #self.p_ =  self.pAS.p()/1000;
        #self.rho_ = self.pAS.rhomass()
    
    cpdef double s_liq(self, string Liq, double T) except *:
        self.T0 = 273.15
        self.P0 = 101.325

        if Liq == b"PAO":
            self.s_l = 1.940 * log(T/self.T0)
        elif Liq == b"PAG":
            try:
                self.s_l = 2.74374E-03*(T-self.T0)+1.08646*log(T/self.T0)
            except:
                a=4
        elif Liq == b"POE":
            self.s_l = 2.30 * log(T/self.T0)
        elif Liq == b"Duratherm_LT":
            self.s_l = (3.4014*(T-298)+1094.3*log(T/298.0))/1000 
        elif Liq == b"Water":
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
    
    cpdef double s_mix(self, string Ref, string Liq, double T, double P, double xL) except *: 
        """
        Entropy of the mixture as a function of temperature [K] and pressure [kPa].  Output in kJ/kg-K
        """
        self.s_g = self.pAS.keyed_output(CoolProp.iSmass)/1000.0
        self.s_m = xL*self.s_liq(Liq,T) + (1-xL)*self.s_gas(Ref,T,P)
        
        if isnan(self.s_m) == True:
            print 's_m is a NaN'
        if isinf(self.s_m) == True:
            print 's_m is Infinite'
        else:
            return self.s_m
        
    cpdef double u_liq(self, string Liq, double T) except *:
        """
        """
        if Liq == b"Duratherm_LT":
            #internal energy [kJ/kg] of Duratherm LT given T in K
            self.u_l = (3.4014/2*pow(T,2)+1094.3*T)/1000  #LT
        elif Liq == b"Zerol60":
            self.u_l = (5.186/2*pow(T,2)+337.116*T)/1000  #Zerol60
            
        elif Liq == b"POE":
            self.u_l = (2.0935/2*T**2 + 1186.7*T)/1000 #POE 150 SUS
            
        else:
            print "Invalid Fluid"
        
        return self.u_l*1000  #J/kg
    
    cpdef double u_mix(self,string Ref, string Liq, double T, double P, double xL) except *:
        """
        input in K, [-] output in J/kg
        """
        
        self.u_l = self.u_liq(Liq,T)
        self.u_g = self.pAS.keyed_output(CoolProp.iUmass)/1000.0
        
        self.u_m = xL*self.u_l + (1-xL)*self.u_g
        
        if isnan(self.u_m) == True:
            print 'u_m is a NaN'
        if isinf(self.u_m) == True:
            print 'u_m is Infinite'
        else:
            return self.u_m  #J/kg
    
    cpdef double h_liq(self, string Liq, double T, double P) except *:  
        self.T0 = 273.15
        self.P0 = 101.325
        self.h = 0
        self.h_l = 0
    
        if Liq == b'PAO':
            self.h_l = 1.940*(T-self.T0)+(P-self.P0)/849
        elif Liq == b'PAG':
            # PAG 0-OB-1020 from Tribology Data Handbook
            self.rho_l = -0.726923*float64(T)+1200.22;
            self.h_l = 2.74374E-03*(float64(T)**2-self.T0**2)/2.0+1.08646*(float64(T)-self.T0)+(float64(P)-self.P0)/self.rho_l;
        elif Liq == b'POE':
            # From Totten, p 261, cp=0.55 cal/g-C --> 2.30 kJ/kg-K
            self.h_l = 2.30*(T-self.T0)+(P-self.P0)/930
        elif Liq == b'Duratherm_LT':
            #the specific enthalpy of Duratherm LT [kJ/kg-k]"
            self.h_l = self.u_liq(Liq,T)/1000 + (P-self.P0)/self.rho_liq(Liq,T)
        elif Liq == b'Water':
            cl_A = 92.053
            cl_B = -0.039953
            cl_C = -0.00021103
            cl_D = 5.3469E-07
            MM_l = 18.0153
            self.h_l = (cl_A*(T-298.15) + cl_B/2.0*(T**2-298.15**2) + cl_C/3.0*(T**3-298.15**3) + cl_D/4.0*(T**4-298.15**4))/MM_l+(P-self.P0)/self.rho_liq(Liq,T)
        else:
            print "Invalid fluid"
            
        return self.h_l*1000
        
    cpdef double h_mix(self, string Ref, string Liq, double T, double P, double xL) except *:
        
        self.h_l = self.h_liq(Liq,T,P)
        self.h_g = self.pAS.keyed_output(CoolProp.iHmass)/1000.0
        
        self.h_m = xL*self.h_l + (1-xL)*self.h_g
    
        if isnan(self.h_m) == True:
            print 'h_m is a NaN'
        if isinf(self.h_m) == True:
            print 'h_m is Infinite'
        else:
            return self.h_m
        
    cpdef double rho_liq(self, string Liq, double T) except *:
        if Liq == b'PAO':
            self.rho_l = 849
        elif Liq == b'PAG':  # PAG 0-OB-1020 from Tribology Data Handbook
            self.rho_l = -0.726923*T+1200.22;
        elif Liq == b'POE':
            self.rho_l = 930
        elif Liq == b'Duratherm_LT':
        #density [kg/m^3] of Duratherm LT given T in K"
            self.rho_l = -0.6793*T + 1012.4 
        elif Liq == b"Water":
            # Water Props from Yaws
            rhol_A=0.3471     
            rhol_B=0.274      
            rhol_n=0.28571    
            rhol_Tc=647.13
            self.rho_l = rhol_A/pow(rhol_B,pow(1-T/rhol_Tc,rhol_n))*1000;
        else:
            print "Invalid fluid"    
        return self.rho_l
    
    cpdef double rho_mix(self, string Ref, string Liq, double T, double P, double xL) except *:
        cdef double rho, v_g, v_l, x_g, S
        
        v_l = 1.0/self.rho_liq(Liq, T)
        v_g = 1.0/self.pAS.rhomass()
        x_g = 1.0 - xL
        x = x_g
        
        if self.model == 'HEM':
            S = 1
        elif self.model == 'Zivi':
            S = (self.v_g/self.v_l)**(0.33333) # Eqn. 4.57 from Chisholm
        elif self.model == 'Fauske':
            S = (self.v_g/self.v_l)**(0.5)     # Eqn. 4.57 from Chisholm
        
        rho = (x + S*(1.0-x))/(x*self.v_g + S*(1.0-x)*self.v_l) #Eq 2.36 from Chisholm
        
        if x>0:
            self.alpha = 1.0/(1.0+(1.0-x)/x*self.rho_g/self.rho_l*self.S)
        else:
            self.alpha = 0.0

        return rho

    cpdef double cp_liq(self, string Liq, double T) except *:
    
        if Liq == b'PAO':
            self.cp_l = 1.940
        elif Liq == b'POE':
            self.cp_l = 2.30
        elif Liq == b'PAG':
            # PAG 0-OB-1020 from Tribology Data Handbook
            # T in K, cp in kJ/kg-K
            self.cp_l = 2.74374E-03*T+1.08646;
        elif Liq == b'Duratherm_LT':
            #specific heat [kJ/kg-K] of Duratherm LT given T in K
            self.cp_l = (3.4014*T + 1094.3)/1000
        elif Liq == b"Water":
            MM_l=18.0153
            cl_A=92.053/MM_l
            cl_B=-0.039953/MM_l
            cl_C=-0.00021103/MM_l
            cl_D=5.3469E-07/MM_l
            self.cp_l = cl_A + cl_B*T + cl_C*T*T + cl_D*T*T*T
        elif Liq == b"ACD100FY":
            # 273 < T [K] < 387
            self.cp_l = 1.304 + 1.035e-3*T+2.801e-6*T**2   #[kJ/kg-K]
        return self.cp_l*1000
            
    cpdef double cp_mix(self,string Ref, string Liq, double T, double P, double xL) except *:
        
        xL = float(xL)
        self.cp_g = self.pAS.cpmass()/1000.0
        self.cp_l = self.cp_liq(Liq,T)
        
        self.cp_m = xL*self.cp_l + (1.0-xL)*self.cp_g
    
        if isnan(self.cp_m) == True:
            print 'cpm is a NaN'
        if isinf(self.cp_m) == True:
            print 'cpm is Infinite'
        else:
            return self.cp_m
        
    cpdef double cv_mix(self, string Ref, string Liq, double T, double P, double xL) except *:
        
        xL = float(xL)
        self.cv_g = self.pAS.cvmass()/1000.0
        self.cv_l = self.cp_liq(Liq,T)
        
        self.cv_m = xL*self.c_l + (1.0-xL)*self.cv_g
    
        if isnan(self.cv_m) == True:
            print 'cvm is a NaN'
        if isinf(self.cv_m) == True:
            print 'cvm is Infinite'
        else:
            return self.cv_m
    
    cpdef double mu_liq(self, string Liq, double T) except *:
        """
        Dynamic Viscosity [Pa-s]
        """
        self.mu_l = 0
        
        if Liq == b'Duratherm_LT': 
            self.mu_l = 8e12*pow(T,-6.001)  #LT
        elif Liq == b'Zerol60': 
            self.mu_l = 1.0*(-0.0001235*T + 0.04808) #Zerol60 Vincent
        elif Liq == b'POE':
            #POE equation only valid from 60 C to 120 C
            #POE 150 SUS
            self.mu_l = 0.000000000517*T**4 - 0.000000795840*T**3 + 0.000460766590*T**2 - 0.118976538068*T + 11.571730524692                   
        elif Liq == b"ACD100FY":
            # 313 < T [K] < 387
            #mu_l = 1.603e+46*pow(T,-19.69) #T [K]   
            self.mu_l = 2.0022e+31*pow(T,-1.2961e+01)
        else:
            print "Invalid fluid"
            
        return self.mu_l
        
    cpdef double mu_mix(self, string Ref, string Liq, double T, double P, double xL) except *:
        self.mu_l = self.mu_liq(Liq,T)
        self.mu_g = self.pAS.keyed_output(CoolProp.iviscosity)
        
        self.mu_m = 1/(xL/self.mu_l + (1-xL)/self.mu_g )
        if isnan(self.mu_m) == True:
            print 'mu_m is a NaN'
        if isinf(self.mu_m) == True:
            print 'mu_m is Infinite'
        else:
            return 1/(xL/self.mu_l + (1-xL)/self.mu_g )
            
    cpdef double VoidFrac(self, string Ref, string Liq, double T, double P, double xL) except *:
        """
        Input in K,kPa , [-] output [-] 
        """
        self.rho_l = self.rho_liq(Liq,T)
        self.rho_g = self.rho_gas(Ref,T,P)
        self.VF = ((1-xL)/self.rho_g)/((1-xL)/self.rho_g + xL/self.rho_l)
        
        if isnan(self.VF) == True:
            print 'VF is a NaN'
        if isinf(self.VF) == True:
            print 'VF is Infinite'
        else:
            return self.VF
    
    cpdef double k_liq(self, string Liq, double T) except *:
    
        if Liq == b'Duratherm_LT': 
            self.k_l = -9e-5*T + 0.1223
        elif Liq == b'Zerol60': 
            self.k_l =  0.1153 #Zerol60 Vincent
        elif Liq == b'POE':  
            self.k_l = 0.138   #POE 150 SUS               
        elif Liq == b"ACD100FY":
            self.k_l = 0.13    #TODO: check value of conductivity
        else:
            print "Invalid fluid"
            
        return self.k_l
        
    cpdef double k_mix(self, string Ref, string Liq, double T, double P, double xL) except *:
        """
        Input in K, kPa , [-] output in kW/m-K
        """
        self.k_l = self.k_liq(Liq,T)
        self.k_g = self.pAS.keyed_output(CoolProp.iconductivity)
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
       
    def dudxL_mix(self,Ref,Liq,T,P,xL):

        delta =.001;
        self.dudxL_m = (self.u_mix(Ref,Liq,T,P,xL+delta) - self.u_mix(Ref,Liq,T,P,xL))/delta
        
        return self.dudxL_m
        
    def T_sp(self,Ref,Liq,s,p,xL,T_guess):   
        """
        Solve for the temperature which gives the same entropy - s [kJ/kg-K], p [kPa], T [K]
        """
    
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
        """
        Solve for the temperature which gives the same enthalpy - h [kJ/kg], p [kPa], T [K]
        """
    
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
    
    def T_crit(self):    
        """
        Return the critical temperature of the refrigerant [K]
        """
        return self.pAS.keyed_output(constants_header.iT_critical)
    
    def p_crit(self):    
        """
        Return the critical pressure of the refrigerant [kPa]
        """
        return self.pAS.keyed_output(constants_header.iP_critical)/1000.0
        
    cpdef double get_Q_m(self) except *:
        """ Get the quality [-] """
        return self.Q
    property Q:
        """ The quality [-] """
        def __get__(self):
            return self.get_Q()

    cpdef double get_rho(self) except *:
        """ Get the density [kg/m^3] """
        return self.rho_mix(self.Ref, self.Liq, self.T, self.p, self.xL)
    property rho:
        """ The density [kg/m^3] """
        def __get__(self):
            return self.get_rho()

    cpdef double get_p(self) except *:
        """ Get the pressure [kPa] """
        return self.P
    property p:
        """ The pressure [kPa] """
        def __get__(self):
            return self.get_p()

    cpdef double get_T(self) except *:
        """ Get the temperature [K] """
        return self.T
    property T:
        """ The temperature [K] """
        def __get__(self):
            return self.get_T()

    cpdef double get_h(self) except *:
        """ Get the specific enthalpy of the mixture [kJ/kg] """
        return self.h_mix(self.Ref, self.Liq, self.T, self.p, self.xL)
    property h:
        """ The specific enthalpy [kJ/kg] """
        def __get__(self):
            return self.get_h()

    cpdef double get_u(self) except *:
        """ Get the specific internal energy [kJ/kg] """
        return self.u_mix(self.Ref, self.Liq, self.T, self.p, self.xL)
    property u:
        """ The internal energy [kJ/kg] """
        def __get__(self):
            return self.get_u()

    cpdef double get_s(self) except *:
        """ Get the specific enthalpy [kJ/kg/K] """
        return self.s_mix(self.Ref, self.Liq, self.T, self.p, self.xL)
    property s:
        """ The specific enthalpy [kJ/kg/K] """
        def __get__(self):
            return self.get_s()

    cpdef double get_cp(self) except *:
        """ Get the specific heat at constant pressure  [kJ/kg/K] """
        return self.cp_mix(self.Ref, self.Liq, self.T, self.p, self.xL)
    property cp:
        """ The specific heat at constant pressure  [kJ/kg/K] """
        def __get__(self):
            return self.get_cp()

    cpdef double get_cv(self) except *:
        """ Get the specific heat at constant volume  [kJ/kg/K] """
        return self.cv_mix(self.Ref, self.Liq, self.T, self.p, self.xL)
    property cv:
        """ The specific heat at constant volume  [kJ/kg/K] """
        def __get__(self):
            return self.get_cv()

    cpdef double get_visc(self) except *:
        """ Get the viscosity, in [Pa-s]"""
        return self.mu_mix(self.Ref, self.Liq, self.T, self.p, self.xL)
    property visc:
        """ The viscosity, in [Pa-s]"""
        def __get__(self):
            return self.get_visc()

    cpdef double get_cond(self) except *:
        """ Get the thermal conductivity, in [kW/m/K]"""
        return self.k_mix(self.Ref, self.Liq, self.T, self.p, self.xL)
    property k:
        """ The thermal conductivity, in [kW/m/K]"""
        def __get__(self):
            return self.get_cond()

    cpdef double get_dpdT(self) except *:
        raise ValueError()
        return self.dpdT_const_V
    property dpdT:
        def __get__(self):
            return self.get_dpdT()
             
    cpdef double get_dudxL(self) except *:
        raise ValueError()
        return self.dudxL_mix
    property dudxL:
        def __get__(self):
            return self.get_dudxL()
            