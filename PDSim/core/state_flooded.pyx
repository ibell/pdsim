from __future__ import division

cimport cython
import CoolProp.CoolProp
import math
from math import floor
from math import log
from math import sqrt
import warnings
from scipy import optimize
from numpy import float64,isnan,isinf,fabs
from scipy.optimize import newton,fsolve

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
        self.p_ = P # Defined in base class
        self.T_ = T # Defined in base class
        self.Liq = Liq # Defined in this class
        self.xL_ = xL # Defined in this class
    
        self.model = model
        
        if b'::' in Ref:
            backend, Ref = Ref.split(b'::', 1)

        self.set_Fluid(Ref, backend)
        self.phase = b''

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
    
    cpdef update(self, dict params):
        """
        Parameters
        params, dictionary
            A dictionary of terms to be updated, with keys equal to single-char inputs to the Props function,
            for instance ``dict(T=298, P = 101.325)`` would be one standard atmosphere
        """

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
        
        self.T_  = params['T']
        self.p_  = params['P']
        self.xL_ = params['xL']
        self.pAS.update(constants_header.PT_INPUTS, self.p_*1000, self.T_)
        #self.rho_ = self.pAS.rhomass()
        
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
    
    cpdef StateFlooded copy2(self):
        """
        Make a copy of this StateFlooded class
        """
        cdef StateFlooded S = StateFlooded(self.Fluid,self.Liq,self.p_,self.T_,self.xL_,self.model)
        #S.phase = self.phase
        return S
    
    
    cpdef double s_liq(self) except *:
        """
        Specific entropy of the flooding medium
        Parameters
        ----------
        T: float
            Temperature [K]
        Liq: string
             Current lubricant implemented: "PAO", "PAG", "POE", "Duratherm LT", "Zerol 60","Water"
             
        Returns
        -------
        s_l: float
            specific entropy [kJ/kg-K]
        """
        
        cdef double T0 = 273.15, P0 = 101.325, s_l, T = self.T_

        if self.Liq == b"PAO":
            s_l = 1.940 * log(T/T0)
        elif self.Liq == b"PAG":
            try:
                s_l = 2.74374E-03*(T-T0)+1.08646*log(T/T0)
            except:
                a=4
        elif self.Liq == b"POE":
            s_l = 2.30 * log(T/T0)
        elif self.Liq == b"Duratherm_LT":
            s_l = (3.4014*(T-298)+1094.3*log(T/298.0))/1000 
        elif self.Liq == b"Zerol60":
            s_l = (5.186*(T-298)+3337.116*log(T/298.0))/1000 
        elif self.Liq == b"Water":
            cl_A=92.053
            cl_B=-0.039953
            cl_C=-0.00021103
            cl_D=5.3469E-07
            MM_l=18.0153
            try:
                s_l = (cl_A*log(T/298.15) + cl_B*(T-298.15) + cl_C/2.0*(T*T-298.15*298.15) + cl_D/3.0*(T*T*T-298.15*298.15*298.15))/MM_l
            except:
                a=1
        else:
            raise ValueError(b"Invalid liquid:" + self.Liq)
    
        return s_l #kJ/kg-K
    
    cpdef double s_mix(self) except *: 
        """
        Entropy of the mixture 
        
        Returns
        ---------
        s_m: float
            specific entropy [kJ/kg-K]
        """
        
        cdef double s_g, s_l
        s_g = self.pAS.keyed_output(constants_header.iSmass)/1000.0
        s_m = self.xL_*self.s_liq() + (1-self.xL_)*s_g
        return s_m
        
    cpdef double u_liq(self) except *:
        """
        Specific internal energy of the flooding medium
        Parameters
        ----------
        T: float
            Temperature [K]
        Liq: string
             Current lubricant implemented: "PAO", "PAG", "POE", "Duratherm LT", "Zerol 60","Water"
             
        Returns
        -------
        u_l: float
            specific internal energy [kJ/kg]
        """
        
        cdef double T = self.T_
        if self.Liq == b"Duratherm_LT":
            #internal energy [kJ/kg] of Duratherm LT given T in K
            u_l = (3.4014/2*T**2+1094.3*T)/1000  #LT
        elif self.Liq == b"Zerol60":
            u_l = (5.186/2*T**2+337.116*T)/1000  #Zerol60
        elif self.Liq == b"POE":
            u_l = (2.0935/2*T**2 + 1186.7*T)/1000 #POE 150 SUS
        else:
            raise ValueError(b"Invalid liquid:" + self.Liq)
        
        return u_l  #kJ/kg
    
    cpdef double u_mix(self) except *:
        """
        Specific internal energy of the mixture 
        
        Returns
        ---------
        u_m: float
            specific internal energy [kJ/kg]
        """
        
        cdef double u_l, u_g, u_m
        u_l = self.u_liq()
        u_g = self.pAS.keyed_output(constants_header.iUmass)/1000.0
        u_m = self.xL_*u_l + (1-self.xL_)*u_g
        return u_m
    
    cpdef double h_liq(self) except *:  
        """
        Specific enthalpy of the flooding medium
        Parameters
        ----------
        T: float
            Temperature [K]
        P: float
            Pressure [kPa]
        Liq: string
             Current lubricant implemented: "PAO", "PAG", "POE", "Duratherm LT", "Zerol 60","Water"
             
        Returns
        -------
        h_l: float
            specific enthalpy[kJ/kg]
        """        
        
        
        cdef double T0 = 273.15, P0 = 101.325, h = 0, h_l = 0, T = self.T_, P = self.p_
        if self.Liq == b'PAO':
            h_l = 1.940*(T-T0)+(P-P0)/849
        elif self.Liq == b'PAG':
            # PAG 0-OB-1020 from Tribology Data Handbook
            rho_l = -0.726923*T+1200.22;
            h_l = 2.74374E-03*(T**2-T0**2)/2.0+1.08646*(T-T0)+(P-P0)/self.rho_liq();
        elif self.Liq == b'POE':
            # From Totten, p 261, cp=0.55 cal/g-C --> 2.30 kJ/kg-K
            h_l = 2.30*(T-T0)+(P-P0)/930
        elif self.Liq == b'Duratherm_LT':
            #the specific enthalpy of Duratherm LT [kJ/kg-k]"
            h_l = self.u_liq()/1000 + (P-P0)/self.rho_liq()
        elif self.Liq == b"Zerol60":
            h_l = self.u_liq()/1000 + (P-P0)/self.rho_liq()
        elif self.Liq == b'Water':
            cl_A = 92.053
            cl_B = -0.039953
            cl_C = -0.00021103
            cl_D = 5.3469E-07
            MM_l = 18.0153
            h_l = (cl_A*(T-298.15) + cl_B/2.0*(T**2-298.15**2) + cl_C/3.0*(T**3-298.15**3) + cl_D/4.0*(T**4-298.15**4))/MM_l+(P-P0)/self.rho_liq()
        else:
            raise ValueError(b"Invalid liquid:" + self.Liq)
            
        return h_l #kJ/kg
        
    cpdef double h_mix(self) except *:
        """
        Specific enthalpy of the mixture 
        
        Returns
        ---------
        h_m: float
            specific internal energy [kJ/kg]
        """
        
        cdef double h_g, h_m
        h_g = self.pAS.keyed_output(constants_header.iHmass)/1000.0
        h_m = self.xL_*self.h_liq() + (1-self.xL_)*h_g
        return h_m

    cpdef double e_mix(self) except *:
       """
       Specific Exergy of the mixture
       
       
       Returns
       -------
       e_mix: float
                Specific exergy [kJ/kg]
       """
       
       cdef double T0_ref = 25+273.15, P0_ref = 101.325, 
       cdef double e_m, h_mix, h_mix_ref, s_mix, s_mix_ref
  
       h_mix = self.h_mix()
       s_mix = self.s_mix()
       
       self.update(dict(T = T0_ref, P = P0_ref, xL = self.xL_))
       
       h_mix_ref = self.h_mix()
       s_mix_ref = self.s_mix()
       
       
       e_m = h_mix - h_mix_ref - T0_ref*(s_mix) # - s_mix_ref
    
       if isnan(e_m) == True:
           print 'e_m is a NaN'
       if isinf(e_m) == True:
           print 'e_m is Infinite'
       else:
           return e_m


    cpdef double rho_liq(self) except *:
        """
        Density of the flooding medium
        Parameters
        ----------
        T: float
            Temperature [K]
        Liq: string
             Current lubricant implemented: "PAO", "PAG", "POE", "Duratherm LT", "Zerol 60","Water"
             
        Returns
        -------
        rho_l: float
            density [kg/m^3]
        """        
        
        cdef double rho_l, T = self.T_
        if self.Liq == b'PAO':
            rho_l = 849
        elif self.Liq == b'PAG':  # PAG 0-OB-1020 from Tribology Data Handbook
            rho_l = -0.726923*T+1200.22;
        elif self.Liq == b'POE':
            rho_l = 930
        elif self.Liq == b'Duratherm_LT':
            #density [kg/m^3] of Duratherm LT given T in K"
            rho_l = -0.6793*T + 1012.4 
        elif self.Liq == b"Zerol60":
            rho_l = -0.667*T + 1050.86
        elif self.Liq == b"Water":
            # Water Props from Yaws
            rhol_A=0.3471     
            rhol_B=0.274      
            rhol_n=0.28571    
            rhol_Tc=647.13
            rho_l = rhol_A/pow(rhol_B,pow(1-T/rhol_Tc,rhol_n))*1000;
        else:
            raise ValueError(b"Invalid liquid:" + self.Liq)    
        return rho_l
    
    cpdef double rho_mix(self) except *:
        """
        Density of the mixture 
        
        Parameters
        ----------
        model: string
                Slip ratio correlation: "HEM", "Zivi","Fauske"
        
        Returns
        ---------
        rho_m: float
            specific internal energy [kJ/kg]
        
        """
        
        cdef double rho, v_g, v_l, x_g, S, x
        
        v_l = 1.0/self.rho_liq()
        v_g = 1.0/self.pAS.rhomass()
        x_g = 1.0 - self.xL_
        x = x_g
        
        if self.model == 'HEM':
            S = 1
        elif self.model == 'Zivi':
            S = (v_g/v_l)**(0.33333) # Eqn. 4.57 from Chisholm
        elif self.model == 'Fauske':
            S = (v_g/v_l)**(0.5)     # Eqn. 4.57 from Chisholm
        
        rho_m = (x + S*(1.0-x))/(x*v_g + S*(1.0-x)*v_l) #Eq 2.36 from Chisholm
        
        if x>0:
            self.alpha = 1.0/(1.0+(1.0-x)/x*v_l/v_g*S)
        else:
            self.alpha = 0.0
        return rho_m

    cpdef double cp_liq(self) except *:
        """
        Specific heat of the flooding medium
        Parameters
        ----------
        T: float
            Temperature [K]
        Liq: string
             Current lubricant implemented: "PAO", "PAG", "POE", "Duratherm LT", "Zerol 60","Water"
             
        Returns
        -------
        cp_l: float
            specific heat [kJ/kg-K]
        """ 
        
        cdef double T = self.T_, cp_l
        if self.Liq == b'PAO':
            cp_l = 1.940
        elif self.Liq == b'POE':
            cp_l = 2.30
        elif self.Liq == b'PAG':
            # PAG 0-OB-1020 from Tribology Data Handbook
            # T in K, cp in kJ/kg-K
            cp_l = 2.74374E-03*T+1.08646;
        elif self.Liq == b'Duratherm_LT':
            #specific heat [kJ/kg-K] of Duratherm LT given T in K
            cp_l = (3.4014*T + 1094.3)/1000
        elif self.Liq == b"Zerol60":
            cp_l = (5.186*T + 337.116)/1000
        elif self.Liq == b"Water":
            one_over_MM_l=1/18.0153
            cl_A=92.053*one_over_MM_l
            cl_B=-0.039953*one_over_MM_l
            cl_C=-0.00021103*one_over_MM_l
            cl_D=5.3469E-07*one_over_MM_l
            cp_l = cl_A + cl_B*T + cl_C*T*T + cl_D*T*T*T
        elif self.Liq == b"ACD100FY":
            # 273 < T [K] < 387
            cp_l = 1.304 + 1.035e-3*T+2.801e-6*T**2   #[kJ/kg-K]
        return cp_l #[kJ/kg]
            
    cpdef double cp_mix(self) except *:
        """
        Specific heat at constant pressure of the mixture 
        Returns
        ---------
        cp_m: float
            specific heat [kJ/kg-K]
        
        """        
        
        cdef double cp_g, cp_l, xL = self.xL_
        cp_g = self.pAS.cpmass()/1000.0
        cp_l = self.cp_liq()
        cp_m = xL*cp_l + (1.0-xL)*cp_g
        return cp_m
        
    cpdef double cv_mix(self) except *:
        """
        Specific heat at constant volume of the mixture 
        Returns
        ---------
        cv_m: float
            specific heat [kJ/kg-K]
        
        """ 

        cdef double cv_g, cv_l, xL = self.xL_
        cv_g = self.pAS.cvmass()/1000.0
        cv_l = self.cp_liq()
        cv_m = xL*cv_l + (1.0-xL)*cv_g
        return cv_m
    
    cpdef double mu_liq(self) except *:
        """
        Dynamic viscosity of the flooding medium
        Parameters
        ----------
        T: float
            Temperature [K]
        Liq: string
             Current lubricant implemented: "PAO", "PAG", "POE", "Duratherm LT", "Zerol 60","Water"
             
        Returns
        -------
        mu_l: float
            dynamic viscosity [Pa-s]
        """
        
        cdef double T = self.T_, mu_l
        
        if self.Liq == b'Duratherm_LT': 
            mu_l = 8e12*T**(-6.001)  #LT
        elif self.Liq == b'Zerol60': 
            mu_l = 1.0*(-0.0001235*T + 0.04808) #Zerol60 Vincent
        elif self.Liq == b'POE':
            #POE equation only valid from 60 C to 120 C
            #POE 150 SUS
            mu_l = 0.000000000517*T**4 - 0.000000795840*T**3 + 0.000460766590*T**2 - 0.118976538068*T + 11.571730524692                   
        elif self.Liq == b"ACD100FY":
            # 313 < T [K] < 387
            #mu_l = 1.603e+46*pow(T,-19.69) #T [K]   
            mu_l = 2.0022e+31*T**(-1.2961e+01)
        else:
            raise ValueError(b"Invalid liquid:" + self.Liq)
        return mu_l
        
    cpdef double mu_mix(self) except *:
        """
        Dynamic viscosity of the mixture 
        Returns
        ---------
        mu_m: float
            dynamic viscosity [Pa-s]
        
        """ 
        
        cdef double mu_l, mu_g, mu_m, xL = self.xL_
        mu_l = self.mu_liq()
        mu_g = self.pAS.keyed_output(constants_header.iviscosity)
        mu_m = 1/(xL/mu_l + (1-xL)/mu_g )
        return mu_m
            
    cpdef double VoidFrac(self) except *:
        """
        Void Fraction of the mixture
        
        Parameters
        ----------
        T: float
            Temperature [K]
        P:float
            Pressure [kPa]
        xL: float
            Liquid mass fraction [-]
        Returns
        -------
        VF: float
            Void fraction [-]
        """
        
        cdef double rho_l, rho_g, xL = self.xL_, VF
        rho_l = self.rho_liq()
        rho_g = self.pAS.rhomass()
        VF = ((1-xL)/rho_g)/((1-xL)/rho_g + xL/rho_l)
        return VF
    
    cpdef double k_liq(self) except *:
        """
        Thermal Conductivity of the flooding medium
        
        Parameters
        ----------
        T: float
            Temperature [K]
        Returns
        -------
        k_l: float
            Thermal conductivity [W/m-K]
        """        
        

        cdef double T = self.T_
        if self.Liq == b'Duratherm_LT': 
            k_l = -9e-5*T + 0.1223
        elif self.Liq == b'Zerol60': 
            k_l =  0.1153 #Zerol60 Vincent
        elif self.Liq == b'POE':  
            k_l = 0.138   #POE 150 SUS               
        elif self.Liq == b"ACD100FY":
            k_l = 0.13    #TODO: check value of conductivity
        else:
            raise ValueError(b"Invalid liquid:" + self.Liq)
            
        return k_l
        
    cpdef double k_mix(self) except *:
        """
        Thermal Conductivity of the mixture
        
        Parameters
        ----------
        T: float
            Temperature [K]
        Returns
        -------
        k_m: float
            Thermal conductivity [kW/m-K]
        """ 
        
        cdef double k_l, k_g, VF, k_m
        k_l = self.k_liq()
        k_g = self.pAS.keyed_output(constants_header.iconductivity)/1000
        VF = self.VoidFrac()
        k_m = (1-VF)*k_l + VF*k_g
        return k_m
    
    cpdef double Pr_mix(self) except *: 
        """
        Prandtl Number of the mixture
        
        Returns
        -------
        Pr_m: float
                Prandtl Number [-]
        """
        
        cdef double Pr_m = (self.cp_mix()*self.mu_mix())/self.k_mix()
        return Pr_m
            
    cpdef double kstar_mix(self) except *:
        
        """
        Ratio of specific heats defined by Hugenroth (2006)
        
        Parameters
        ----------
        cp_g: float
                Specific heat at constant pressure of the refrigerant [kJ/kg-K]
        cv_g: float
                Specific heat at constant volume of the refrigerant [kJ/kg-K] 
        cp_liq: float
                Specific heat of the flooding medium [kJ/kg-K]
        xL: float
            Liquid mass fraction [-]
        
        Returns
        -------
        kstar_m: float
                    Ratio of specific heats of the mixture [-]
        """
        
        
        cdef double xL = self.xL_, kstar_m, cv_g, cp_g
        cp_g = self.pAS.cpmass()/1000.0
        cv_g = self.pAS.cvmass()/1000.0
        kstar_m =((1-xL)*cp_g + xL*self.cp_liq())/((1-xL)*cv_g + xL*self.cp_liq())
        return kstar_m
    

    cpdef double dpdT_const_V(self) except *:
        """
        Derivative of pressure over the temperature at constant specific
        volume of the mixture
        
        Parameters
        ----------
        P: float
            Pressure [kPa]
        v: float
            Specific volume [m^3/kg]
        T: float
            Temperature [K]
        
        Returns
        -------
        dPdT_v: float
                Derivative [kPa/K]
        
        
        """
    
        
        cdef double f,P1,P2,v1,T,delta,eps,change,x1,x2,x3,y1,y2
        cdef int iter
        #Ref = self.Fluid, Liq = self.Liq,, xL = self.xL, T = self.T_, 
        P1 = self.p_
        v1 = 1.0/self.rho_mix()

        delta=1.0e-5
        eps = 1e-6
        change = 999
        iter = 1
        
        while (iter<=3 or fabs(f)>eps) and iter < 100:
            
            if iter ==1:
                x1 = P1
                P2 = x1
            if iter == 2:
                x2 = P1 + 0.001
                P2 = x2
            if iter > 2:
                P2 = x2
            
            T = self.T_    
            self.update(dict(T = T +delta, P = P2, xL = self.xL_))
            f = (1.0/self.rho_mix()) - v1
            
            if iter==1:
                y1 = f
            if iter > 1:
                y2 = f
                x3 = x2-y2/(y2-y1)*(x2-x1)
                change = fabs(y2/(y2-y1)*(x2-x1))
                y1=y2
                x1=x2
                x2=x3
            
            iter = iter+1
            
            if iter > 50:
                print 'dPdT_const_v not converging'
            
        
        if isnan((P2-P1)/delta) == True:
            print 'dpdt_v is a NaN'
        if isinf((P2-P1)/delta) == True:
            print 'dpdt_v is Infinite'
        else:
            return (P2-P1)/delta
       


       
    cpdef double dudxL_mix(self) except *:
        """
        Derivative of the internal energy of the mixture over the liquid mass fraction
        
        Parameters
        ----------
        P: float
            Pressure [kPa]
        T: float
            Temperature [K]
        um: float
            Specific internal energy of the mixture [kJ/kg]
        xL: float
            Liquid mass fraction [-]
        Returns
        -------
        dudxL_m: float
                Derivative [kJ/kg]
        
        """

        cdef double u_m1,u_m2, delta, xL
        delta = 0.001
        u_m1 = self.u_mix()
        
        xL = self.xL_
        self.update(dict(T = self.T_, P = self.p_, xL = xL + delta ))
        u_m2 = self.u_mix()
        
        return (u_m2 - u_m1)/delta
        

    cpdef double cK_e(self) except*:
        """
        Effective slip ratio of the mixture
        
        Equation taken from page 43, equation 4.51 from Chisholm (1983) for
        liquid entrainment in gas. Value of psi 0.4 is recommended from text
        
        Parameters
        ----------
        x: float
            Mass dryness fraction [-]
        w: float
            entrainment [-]
        
        Returns
        -------
        Ke: float
                Effective slip ratio [-]
        """        
        
        
        cdef double K,flag,x,v_g,v_l,w

        x = 1 - self.xL_
        v_l = 1.0/self.rho_liq()
        v_g = 1.0/self.pAS.rhomass()
        flag = 0
        w = 0
        
        if (x ==0 or x ==1):
            """
            Homogenous flow
            """
            return 1
    
        if (flag >0.9 and flag <1.1):
            """
            Entrainment slip ratio
            From Chisholm (1983) w=0.4
            """            
            K= w +(1.0 -w)* sqrt(( v_g/v_l + w*(1-x)/x) /(1+ w*(1-x)/x))
    
        if (flag >1.9 and flag <2.1):
            """
            Chisholms sqrt(vh/vl)
            """
            K= sqrt(1.0+ x*( v_g /v_l -1.0))
    
        if (flag >2.9 and flag <3.1):
    
            K= pow(v_g/v_l ,0.25*.28) 

        return K
    
    cpdef double cv_e(self) except*:
        """
        Effective specific volume of the mixture
        Returns
        -------
        ve: float
                effective specific volume [-]
        """           

        cdef double ve,flag,x,K_e,K_c,v_l,v_g,w
        
        K_e = self.cK_e()
        x = 1 - self.xL_
        v_l = 1.0/self.rho_liq()
        v_g = 1.0/self.pAS.rhomass()
        flag = 0
        w = 0
        
        if (flag >0.9 and flag <1.1):
            """
            // using 5.48 and 5.49 from Chisholm
            // if psi=0, separated flow results ( flag ==2 )
            // if psi=1, homogeneous flow results ( flag == 5)
            // So basically this form is general and captures all possibilities , sep , hom
            , or entrained
            // should use this one
            //     (    [(1 -psi)^2] ) -1
            // Kc= ( w+ [----------] )
            //     (    [ K_e -psi ] )
            """
            Kc =1.0/( w +((1.0 - w) *(1.0 -w))/( K_e -w))
            ve =(x* v_g + K_e *(1.0 -x)* v_l )*(x +(1.0 - x)/Kc)
    
        if (flag >1.9 and flag <2.1):
            """
            // Equation 5.13 from Chisholm for separated flow
            """
            ve =(x* v_g + K_e *(1.0 -x)* v_l )*(x +(1.0 - x)/ K_e )
    
        if (flag >2.9 and flag <3.1):
            """
            // Equation 2.48 from Chisholm
            """
            ve =(1+ w*(1 -x)/x*v_l/ v_g ) /(1+ w*(1 -x)/x)* v_g;
    
        if (flag >3.9 and flag <4.1):
            """
            // using 15 from Morris
            // If going to use this formula , must change pow () to multiply
            """
            ve =(x* v_g + K_e *(1.0 -x)* v_l )*(x +(1.0 - x)/ K_e *(1+ (K_e -1)**2.0 /( (v_g /v_l)**0.5-1)))
    
        if (flag >4.9 and flag <5.1):
            """
            // homogenous
            """
            ve =(x* v_g +(1.0 - x)* v_l )
    
        return ve    
    
    
    
    
    
    # cpdef T_sp(self,s,p,xL,T_guess):   
    #     """
    #     Solve for the temperature which gives the same entropy - s [kJ/kg-K], p [kPa], T [K]
    #     """
    # 
    #     """
    #     When fsolve() is called, it will use a single element
    #     ndarray but we only want a double to be passed along,
    #     so take the first element of the 1-element array, which
    #     is a 64-bit float
    #     """ 
    #     def obj(T):
    #         self.update(dict(T = T, xL = self.xL_, P = self.p_))
    #         self.s_mix() - s
    #     T = optimize.fsolve(obj, T_guess)
    #     return T
    
    # cpdef T_hp(self,Ref,Liq,h,p,xL,T_guess):   
    #     """
    #     Solve for the temperature which gives the same enthalpy - h [kJ/kg], p [kPa], T [K]
    #     """
    # 
    #     global h_mix
    #     
    #     """
    #     When fsolve() is called, it will use a single element
    #     ndarray but we only want a double to be passed along,
    #     so take the first element of the 1-element array, which
    #     is a 64-bit float
    #     """
    #     def obj(T):
    #         self.update(dict(T = T, xL = self.xL_, P = self.p_))
    #         self.h_mix() - h
    #         
    #     T = optimize.fsolve(obj, T_guess)
    #     return T
    
    cpdef double T_crit(self) except *:    
        """
        Critical temperature of the refrigerant 
        
        Returns
        -------
        Tcrit: float
                Temperature [K]
        """
        return self.pAS.keyed_output(constants_header.iT_critical)
    
    cpdef double p_crit(self) except *:    
        """
        Critical pressure of the refrigerant 
        
        Returns
        -------
        pcrit: float
                Pressure [kPa]
        """
        return self.pAS.keyed_output(constants_header.iP_critical)/1000.0
    ##
        
    cpdef double get_p(self) except *:
        """ Get the pressure [kPa] """
        return self.p_
    property p:
        """ The pressure [kPa] """
        def __get__(self):
            return self.get_p()

    cpdef double get_T(self) except *:
        """ Get the temperature [K] """
        return self.T_
    property T:
        """ The temperature [K] """
        def __get__(self):
            return self.get_T()
            
    cpdef double get_Q(self) except *:
        """ Get the quality [-] """
        return -1
    property Q:
        """ The quality [-] """
        def __get__(self):
            return self.get_Q()

    cpdef double get_rho(self) except *:
        """ Get the density [kg/m^3] """
        return self.rho_mix()
    property rho:
        """ The density [kg/m^3] """
        def __get__(self):
            return self.get_rho()

    cpdef double get_h(self) except *:
        """ Get the specific enthalpy of the mixture [kJ/kg] """
        return self.h_mix()
    property h:
        """ The specific enthalpy [kJ/kg] """
        def __get__(self):
            return self.get_h()

    cpdef double get_u(self) except *:
        """ Get the specific internal energy [kJ/kg] """
        return self.u_mix()
    property u:
        """ The internal energy [kJ/kg] """
        def __get__(self):
            return self.get_u()

    cpdef double get_s(self) except *:
        """ Get the specific enthalpy [kJ/kg/K] """
        return self.s_mix()
    property s:
        """ The specific enthalpy [kJ/kg/K] """
        def __get__(self):
            return self.get_s()

    
    cpdef double get_e(self) except *:
        """ Get specific exergy [kJ/kg] """
        return self.e_mix()
    property e:
        """ The specific exergy [kJ/kg] """
        def __get__(self):
            return self.get_e()


    cpdef double get_cp(self) except *:
        """ Get the specific heat at constant pressure  [kJ/kg/K] """
        return self.cp_mix()
    property cp:
        """ The specific heat at constant pressure  [kJ/kg/K] """
        def __get__(self):
            return self.get_cp()

    cpdef double get_cv(self) except *:
        """ Get the specific heat at constant volume  [kJ/kg/K] """
        return self.cv_mix()
    property cv:
        """ The specific heat at constant volume  [kJ/kg/K] """
        def __get__(self):
            return self.get_cv()

    cpdef double get_visc(self) except *:
        """ Get the viscosity, in [Pa-s]"""
        return self.mu_mix()
    property visc:
        """ The viscosity, in [Pa-s]"""
        def __get__(self):
            return self.get_visc()

    cpdef double get_cond(self) except *:
        """ Get the thermal conductivity, in [kW/m/K]"""
        return self.k_mix()
    property k:
        """ The thermal conductivity, in [kW/m/K]"""
        def __get__(self):
            return self.get_cond()

    cpdef double get_kstar(self) except *:
        """Get ratio of specific heats of the mixture, in [-] """
        return self.kstar_mix()
    property kstar:
        def __get__(self):
            return self.get_kstar()
    
    
    cpdef double get_dpdT(self) except *:
        """ Get dpdT_const_V, in [ ]"""
        return self.dpdT_const_V()
    property dpdT:
        def __get__(self):
            return self.get_dpdT()
             
    cpdef double get_dudxL(self) except *:
        return self.dudxL_mix()
    property dudxL:
        def __get__(self):
            return self.get_dudxL()


    cpdef double get_cKe(self) except *:
        return self.cK_e()
    property cKe:
        def __get__(self):
            return self.get_cKe()
    
    cpdef double get_cve(self) except *:
        return self.cv_e()
    property cve:
        def __get__(self):
            return self.get_cve()
