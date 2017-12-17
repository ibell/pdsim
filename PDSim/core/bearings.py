from __future__ import print_function, absolute_import

from math import pi,sqrt,log
import numpy as np
from ._bearings import calculate_epsilon_short, calculate_epsilon_long

class struct: pass

def thrust_bearing(**kwargs):
    """
    Analysis for the journal bearing
    
    We only take keyword arguments to avoid problems with order of parameters
    Parameters must be specified by name!
    
    Parameters
    ----------
    mu : float
        Friction coefficient [-]
    V : float
        Contact velocity [m/s]        
    N : float
        Thrust force [N]
    
    Notes
    -----
    To derive the speed of contact of the orbiting scroll
    
    .. math::
    
        x = {r_o}\cos \left( {{\phi _{ie}} - \\frac{\pi }{2} - \\theta } \\right)

    .. math::
    
        y = {r_o}\sin \left( {{\phi _{ie}} - \\frac{\pi }{2} - \\theta } \\right)
    
    .. math::
    
        \\frac{{dx}}{{d\\theta }} = {r_o}\sin \left( {{\phi _{ie}} - \\frac{\pi }{2} - \\theta } \\right)\omega 
    
    .. math::
    
        \\frac{{dy}}{{d\\theta }} =  - {r_o}\cos \left( {{\phi _{ie}} - \\frac{\pi }{2} - \\theta } \\right)\omega

    .. math::
    
        \left| v \\right| = \sqrt {{{\left( {\\frac{{dx}}{{d\\theta }}} \\right)}^2} + {{\left( {\\frac{{dy}}{{d\\theta }}} \\right)}^2}}  = {r_o}\omega
        
    But it is quite possible that you do not have hydro-dynamic lubrication, and 
    as a result you can get asperity-asperity contact and much higher friction
    coefficients - see for example Kobayashi et al. "Experimental Study on Journal Bearing
    Characteristics in Reciprocating Compressors for HFC-134a" Purdue Compressor
    Conferences 1998.  http://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=2257&context=icec
    
    Can be off by a factor of as much as 10 times at low
    Sommerfeld number
     
    """
    #Friction coefficient
    mu = kwargs.pop('mu', None)
    V = kwargs.pop('V', None)
    N = kwargs.pop('N', None)
    
    Wdot_loss_thrust = mu * V * N
    return dict(mu = mu,
                V = V,
                Wdot_loss = Wdot_loss_thrust
                )

def journal_bearing(**kwargs):
    
    """
    The necessary calculations for a journal bearing
    
    We only take keyword arguments to avoid problems with order of parameters
    
    Parameters must be specified by name!
    
    Parameters
    ----------
    r_b : float
        Radius of journal [m]
    L : float
        Length of journal [m]
    omega : float
        Rotational speed [rad/s]
    eta_0 : float, or 1d numpy array
        Viscosity of lubricant [Pa-s]
    W : float
        Applied load [N]
    c : float
        Bearing clearance [m]
    
    Returns
    -------
    output_dict : dictionary
        Dictionary of output terms, 
        
    Notes
    -----
    
    Short-bearing analysis is used.  Here we are implementing the equations from 
    the Ph.D. dissertation of Jay Kim (Purdue, 2005).  Knowing :math:`w_r` and 
    the dimensions we can then obtain :math:`\epsilon` using a 1D secant solver
    with the initial guess of :math:`\epsilon=0.5`.  Then we obtain 
    :math:`F_{shear}` for the given value of :math:`\epsilon`. 
    
    .. math::
    
        W_r = \\frac{w_r}{\\eta_0\\omega r_b L}\\left(\\frac{c}{L}\\right)^2 = \\frac{\\epsilon}{4(1-\\epsilon^2)^2}[16\\epsilon^2+\\pi^2(1-\\epsilon^2)]^{1/2}

    .. math::
    
        \mathbf{F}_{shear} = \\frac{F_{shear}}{\\eta_0\\omega r_b L}\\left(\\frac{c}{r_b}\\right) = \\frac{2\\pi}{(1-\\epsilon^2)^{1/2}}\\left[\\left(\\frac{b}{r_b}\\right)^2\\frac{\\epsilon^2}{16(1-\\epsilon^2)}+1\\right]
        
    .. math::
    
        \\mu = \\frac{F_{shear}}{w_r}

    """

    L = kwargs['L']
    r_b = kwargs['r_b']
    W = kwargs['W']
    eta_0 = kwargs['eta_0']
    omega = kwargs['omega']
    c = kwargs['c']
       
    # Short bearing analysis
    Wr_short = W/(eta_0*omega*r_b*L)*(c/L)**2
    
    # Long bearing analysis
    Wr_long = W/(eta_0*omega*r_b*L)*(c/r_b)**2
    
    #Select the first guess value based on the input to the function
    if Wr_long < 5:
        epsilon0_long = 0.01
    else:
        epsilon0_long = 0.8
        
    if Wr_short < 0.8:
        epsilon0_short = 0.01
    else:
        epsilon0_short = 0.8
    
    # Sommerfeld number
    S = eta_0*omega*r_b*L/(pi*W)*(r_b/c)**2
    
    if L/(2*r_b) < 0.5: #Short bearing analysis
        epsilon = calculate_epsilon_short(np.log(Wr_short), epsilon0_short)
        Fshear = 2*pi/np.sqrt(1-epsilon**2)*((L/r_b)**2*epsilon**2/(16*(1-epsilon**2)) +1)*eta_0*omega*r_b*L*(r_b/c)
    
    elif L/(2*r_b) > 2.0: #Infinitely-long analysis
        epsilon = calculate_epsilon_long(np.log(Wr_long), epsilon0_long)
        Fshear = pi/np.sqrt(1-epsilon**2)*(5*epsilon**2+4)/(epsilon**2+2)*eta_0*omega*r_b*L*(r_b/c)
        
    else: #Weighted average of the two to calculate the shear force
        epsilon = calculate_epsilon_long(np.log(Wr_long),epsilon0_long)
        if not (epsilon >= 0 and epsilon <= 1):
            #print 'Wr_long,epsilon,W',Wr_long,epsilon0_long,W
            if epsilon < 0:
                epsilon = 0
            elif epsilon > 1:
                raise ValueError('epsilon [{epsilon:g}] is greater than 1 for W of [{W:g}] N.  Is your journal bearing too small?'.format(epsilon = epsilon, W = W))
        Fshear_long = pi/np.sqrt(1-epsilon**2)*(5*epsilon**2+4)/(epsilon**2+2)*eta_0*omega*r_b*L*(r_b/c)
        epsilon = calculate_epsilon_short(np.log(Wr_short),epsilon0_short)
        assert epsilon >= 0 and epsilon <= 1
        Fshear_short = 2*pi/np.sqrt(1-epsilon**2)*((L/r_b)**2*epsilon**2/(16*(1-epsilon**2)) +1)*eta_0*omega*r_b*L*(r_b/c)
        Fshear = (Fshear_long-Fshear_short)/1.5*(L/(2*r_b)-0.5) + Fshear_short
    
    hmin_over_c = 1-epsilon
    
    h_min = hmin_over_c*c
    
    f = Fshear/W
    
    #Frictional losses [W]
    Wdot_loss = omega*r_b*f*W
    
    return dict(
                S = S,
                f = f,
                c = c,
                Load = W,
                r_b = r_b,
                h_min = h_min,
                epsilon = epsilon,
                omega = omega,
                eta_0 = eta_0,
                Wdot_loss = Wdot_loss
                )
         
def journal_bearing_design(**kwargs):
    
    """
    The necessary calculations for a journal bearing
    
    We only take keyword arguments to avoid problems with order of parameters
    
    Parameters must be specified by name!
    
    Parameters
    ----------
    r_b : float
        Radius of journal [m]
    L : float
        Length of journal [m]
    omega : float
        Rotational speed [rad/s]
    eta_0 : float
        Viscosity of lubricant [Pa-s]
    W : float
        Applied load [N]
    c : float
        Bearing clearance [m]
    design : string or float
        Either one of ``'friction'`` (design for minimum friction) 
        ``'load'`` (design for maximum load) or a floating point 
        value in the range [0,1] that weights the minimum friction
        and maximum load parameters and 0 gives the friction solution.
        Hamrock recommends a value of 0.5 for general applications
        
    Notes
    -----
    One of ``design`` or ``c`` must be provided
    
    Based on the method presented by
    
        Raimondi, A. A., and Boyd, J. (1958) : A Solution for the Finite 
        Journal Bearing and Its Application to Analysis and Design-I, -II, 
        and -III. ASLE Trans., vol. 1, no. I, I- pp. 159-174; II- pp. 175-193; 
        III- pp. 194-209.
    
    And further analysis presented in Hamrock:
    
        In Fig. 11.2 a recommended operating eccentricity ratio, or minimum film
        thickness, is indicated as well as a preferred operating area. The left boundary of
        the shaded zone defines the optimum eccentricity ratio for a minimum coefficient
        of friction, and the right boundary the optimum eccentricity ratio for maximum
        load. The recommended operating eccentricity for general application is midway
        between these two boundaries.

    """

    L = kwargs['L']
    r_b = kwargs['r_b']
    W = kwargs['W']
    eta_0 = kwargs['eta_0']
    omega = kwargs['omega']
    c = kwargs.get('c',None)
    design = kwargs.get('design',None)
        
    L_over_D = L/(2*r_b)
    
    #First based on L/D and the design type, select the hm_over_c
    if design == 'friction':
        hm_over_c = 0.032*L_over_D**2 + 0.32*L_over_D - 0.052
    elif design == 'load':
        hm_over_c = 0.188272*log(L_over_D) + 0.541167
    elif isinstance(design,float):
        hm_over_c_friction = 0.032*L_over_D**2 + 0.32*L_over_D - 0.052
        hm_over_c_load = 0.188272*log(L_over_D) + 0.541167
        hm_over_c = design*hm_over_c_load+(1-design)*hm_over_c_friction
    else:
        raise ValueError('invalid value for parameter "design"')
    
    #Interpolate in the data from to obtain S and rb/c*f
    _hm_c = [0.9,0.8,0.6,0.4,0.2,0.1,0.03]
    _L_D = [1,0.5,0.25]
    _S = [[1.33,0.631,0.264,0.121,0.0446,0.0188,0.00474],[4.31,2.03,0.779,0.319,0.0923,0.0313,0.00609],[16.2,7.57,2.83,1.07,0.261,0.0736,0.0101]]
    _rb_c_f = [[26.4, 12.8, 5.79, 3.22, 1.7, 1.05, 0.514],[85.6, 40.9, 17, 8.1, 3.26, 1.6, 0.61], [322, 153, 61.1, 26.7, 8.8, 3.5, 0.922]]

    from scipy.interpolate import interp2d
    f = interp2d(x=_hm_c,
                 y=_L_D,
                 z=_S,
                 kind = 'linear')
    S = f(hm_over_c,L_over_D)[0]
    f = interp2d(x=_hm_c,
                 y=_L_D,
                 z=_rb_c_f,
                 kind = 'linear')
    rb_c_f = f(hm_over_c,L_over_D)[0]
    
    #r_b/c from bearing number
    rb_c = sqrt((S*pi*W)/(eta_0*omega*r_b*L))
    
    #Friction factor after all that [-]
    f = rb_c_f/rb_c 
    
    h_min = hm_over_c*(r_b/rb_c)
    
    c = r_b/rb_c
    
    Wdot_loss = omega*r_b*f*W
    
    return dict(
                S = S,
                f = f,
                c = c,
                Load = W,
                r_b = r_b,
                h_min = h_min,
                epsilon = 1-hm_over_c,
                omega = omega,
                eta_0 = eta_0,
                Wdot_loss = Wdot_loss
                )
            
if __name__ == '__main__':
    d = journal_bearing(r_b = 0.043/2.0,
                        L = 0.057,
                        c = 5e-5,
                        W = 12950,
                        eta_0 = 0.0086,
                        omega = 301.9402420413402)
    
    print(d['epsilon'])
    
    print('The following is the output from example 11.1 in Hamrock')
    print(journal_bearing_design(r_b = 0.02, 
                        L = 0.04, 
                        design = 'load', 
                        W = 2200, 
                        eta_0 = 0.17, 
                        omega = 3600/60.0*2*pi
                        ))
    
