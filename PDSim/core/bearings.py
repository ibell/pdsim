import quantities as pq
from math import pi,sqrt,log

def thrust_bearing(**kwargs):
    """
    
    Notes
    -----
    To derive the speed of contact of the orbiting scroll
    
    .. math::
    
        x = {r_o}\cos \left( {{\phi _{ie}} - \frac{\pi }{2} - \theta } \right)\\

    .. math::
    
        y = {r_o}\sin \left( {{\phi _{ie}} - \frac{\pi }{2} - \theta } \right)\\
    
    .. math::
    
        \frac{{dx}}{{d\theta }} = {r_o}\sin \left( {{\phi _{ie}} - \frac{\pi }{2} - \theta } \right)\omega \\
    
    .. math::
    
        \frac{{dy}}{{d\theta }} =  - {r_o}\cos \left( {{\phi _{ie}} - \frac{\pi }{2} - \theta } \right)\omega \\

    .. math::
    
        \left| v \right| = \sqrt {{{\left( {\frac{{dx}}{{d\theta }}} \right)}^2} + {{\left( {\frac{{dy}}{{d\theta }}} \right)}^2}}  = {r_o}\omega 
    """
    
    pass

def journal_bearing(**kwargs):
    
    """
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
    design : string or float
        Either one of ``'friction'`` (design for minimum friction) 
        ``'load'`` (design for maximum load) or a floating point 
        value in the range [0,1] that weights the minimum friction
        and maximum load parameters and 0 gives the friction solution.
        Hancock recommends a value of 0.5 for general applications
    
    Notes
    -----
    Based on the method presented by
    
        Raimondi, A. A., and Boyd, J. (1958) : A Solution for the Finite Journal Bearing
and Its Application to Analysis and Design-I, -II, and -III. ASLE Trans.,
vol. 1, no. I, I- pp. 159-174; II- pp. 175-193; III- pp. 194-209.
    
    And further analysis presented in 
    From Hancock
    
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
    design = kwargs.get('design','friction')
    
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
    
    #Friction factor after all that
    f = rb_c_f/rb_c
    
    Wdot_loss = omega*r_b*f*W
    
    return dict(
                S = S,
                f = f,
                c = r_b/rb_c,
                Wdot_loss = Wdot_loss
                )
                
if __name__ == '__main__':

    print 'The following is the output from example 11.1 in Hamrock'
    print journal_bearing(r_b = 0.02, 
                        L = 0.04, 
                        design = 'load', 
                        W = 2200, 
                        eta_0 = 0.17, 
                        omega = 3600/60.0*2*pi
                        )
    