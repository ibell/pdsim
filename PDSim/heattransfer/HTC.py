from __future__ import division, print_function
from math import pi, cos, sin, log
import os, sys

#Math Package
from scipy import optimize
from math import ceil
from scipy.special import erf
from PDSim.misc.scipylike import trapz
from scipy.integrate import quad,simps
from scipy.interpolate import interp1d,splrep,splev
import pylab
import numpy as np
# Coolprop
from CoolProp import State,AbstractState
from CoolProp import State
from CoolProp import CoolProp as CP


def HTC(HTCat,T_wall,T_inf,P_film,Ref,L,D_pipe=None,PlateNum=None):

    """
    Nusselt number for different heat transfer categories;
    find heat transfer coefficient based on Nusselt number;
    characteristic length, A/P;

    Return heat transfer rate by radiation
    
    Parameters
    ----------
    A [m^2] : heat transfer area 
    epsilon : surface emissivity
    T_wall [K]: surface temperature
    P_film [kPa]: pressure at the film
    Tinf [K]: surrounding temperature
    sigma [W/(m^2 K)]: Stefan-Boltzmann constant 5.670*10e-8
    ----------
    Return
    h [kW/m^2 K]: Heat transfer coefficient
    """
    g = 9.81 # gravity acceleration [m/s^2]

    # Film temperature, used to calculate thermal propertiers
    T_film = (T_wall + T_inf)/2 # [K]
    
    # thermal expansion coefficient, assumed ideal gas
    beta = 1/T_film # [1/K]

    # Transport properties calcualtion film
    StateFilm = State.State(Ref,dict(T=T_film, P=P_film))   # use the film temperature to find the outer convective coefficient
    
    Pr_film = StateFilm.Prandtl #[-]
    rho_film = StateFilm.rho #[kg/m3]
    k_film = StateFilm.k #[kW/m-K] conductivity
    mu_film = StateFilm.visc #[Pa-s]
    nu_film = mu_film/rho_film  # [m^2/s]
    cp_film = StateFilm.cp # [kJ/kg/K]

    if T_wall<T_inf:
        # Grashof number, beta-thermal expansion coefficient
        Gr = g*beta*(T_inf-T_wall)*L**3/nu_film**2 
    else:
        # Grashof number, beta-thermal expansion coefficient
        Gr = g*beta*(T_wall-T_inf)*L**3/nu_film**2 
    
    # Rayleigh number
    RaL = Gr*Pr_film # [-]
    
    if HTCat == 'horizontal_pipe':
        RaD =  RaL*D_pipe**3/L**3

    if RaL < 1e-3:
        Nu = 0.0
    else:
            if HTCat == 'vertical_plate':
                if RaL > 1e9:
                    Nu = (0.825 + 0.387*RaL**(1/6) / ( 1 + (0.492/Pr_film)**(9/16) )**(8/27))**2
                else:
                    Nu = 0.68 + 0.670*RaL**(1/4) / ( 1 + (0.492/Pr_film)**(9/16) )**(4/9)
                
            elif HTCat == 'horizontal_plate':
                if PlateNum == 'upper_surface': 
                    if T_wall > T_inf: # hot plate
                        if RaL >= 1e4 and RaL <= 1e7:
                            Nu = 0.54*RaL**(1/4)
                        elif RaL >= 1e7 and RaL <= 1e11:
                            Nu = 0.15*RaL**(1/3)         
                        else:
                            Nu = 0.71*RaL**(1/4)    # not quite sure about this equation
                    else: # cold plate
                        if RaL >= 1e5 and RaL <= 1e10:
                            Nu = 0.27*RaL**(1/4)
                        else:
                            Nu = 0.71*RaL**(1/4) # not quite sure about this equation, M. AL-ARABI and M. K. EL-RIEDYt;

                elif PlateNum == 'lower_surface':
                    if T_wall > T_inf: # hot plate
                        if RaL >= 1e5 and RaL <= 1e10:
                            Nu = 0.27*RaL**(1/4)
                        else:
                            Nu = 0.25*RaL**(1/4) # not quite sure about this equation, M. AL-ARABI and M. K. EL-RIEDYt;
                        
                    else: # cold plate
                        if RaL >= 1e4 and RaL <= 1e7:
                            Nu = 0.54*RaL**(1/4)
                        elif RaL >= 1e7 and RaL <= 1e11:
                            Nu = 0.15*RaL**(1/3)
                        else:
                            Nu = 0.25*RaL**(1/4) # not quite sure about this equation, M. AL-ARABI and M. K. EL-RIEDYt;
                else:
                    raise Exception('PlateNum must be either upper_surface or lower_surface')
                
            elif HTCat == 'horizontal_pipe':
                if RaD <= 1e12:
                    # Churchill and Chu, 1975, RaL->RaD    
                    Nu = (0.60+0.387*RaD**(1/6)/(1 + (0.559/Pr_film)**(9/16))**(8/27))**2    

                else: # Kuehn and Goldstein, 1976.
                    temp= (( 0.518*(RaD**0.25)*(1+(0.559/Pr_film)**0.6)**(-5/12) )**15 + (0.1*RaD**(1/3))**15)**(1/15)
                    Nu = 2/(log(1 + 2/temp))
                    
            elif HTCat == 'vertical_pipe':
                if (D_pipe/L) < 35/Gr**(1/4):
                    F = 1/3*((L/D_pipe)/(1/Gr))**(1/4)+1
                else:
                    F = 1.0
                Nu = F*(0.825 + 0.387*RaL**(1/6)/(1+(0.492/Pr_film)**(9/16))**(8/27))**2

            else:
                raise Exception('not implemented')

    # convective heat transfer coefficient
    if HTCat == 'horizontal_pipe':
        h = Nu*k_film/D_pipe 
    else:
        h = Nu*k_film/L

    return h