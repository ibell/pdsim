#cython: embedsignature=True

from __future__ import division

import cython
cimport cython

from libc.math cimport exp, log, M_PI as pi, M_E as e, sqrt

@cython.cdivision(True)
cpdef double p_pstar(double Ma, double gamma):
    """
    Return ``p/pstar``, given by
    
    .. math::
    
        \frac{p}{p^*} = \frac{1}{Ma}\left( \frac{(gamma+1)}{2} \right)^{1/2}*\left(1+\frac{(gamma-1)}{2}Ma^2\right)^{-1/2}
    """
    return 1/Ma*((gamma+1)/2)**0.5*(1+(gamma-1)/2*Ma**2)**(-0.5)

@cython.cdivision(True)
cdef double NondimLength(double Ma, double gamma):
    """
    Returns the logarithm of the non-dimensional length for the Fanno analysis
    Logs are used because of the extreme shape of the objective function
    """
    return log(1/gamma*(1-Ma**2)/Ma**2+(gamma+1)/(2*gamma)*log((1+gamma)*Ma**2/2/(1+(gamma-1)/2*Ma**2)))

@cython.cdivision(True)
cdef double _Ma(double logLparam, double gamma, double Ma0):
    """
    Internal secant method used to calculate Ma for given logarithm of length parameter
    Logs are used because of the extreme shape of the objective function
    """

    cdef double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-10,change=999,f=999,T=300, Ma
    cdef int iter = 1

    while ((iter<=3 or abs(f)>eps) and iter<100):
        if (iter==1):
            x1 = Ma0
            Ma = x1
        if (iter==2):
            x2 = Ma0+0.01
            Ma = x2
        if (iter>2):
            Ma=x2
            
        f = NondimLength(Ma, gamma)-logLparam;
        
        if (iter==1):
            y1=f
        if (iter>1):
            y2=f
            x3=x2-y2/(y2-y1)*(x2-x1)
            change=abs(y2/(y2-y1)*(x2-x1))
            y1=y2; x1=x2; x2=x3
        iter += 1
    return Ma

cpdef double Fanno_Ma_nondimLength(double Lparam, double gamma):
    """
    Fanno flow analysis inverted calculation to get mach number from nondim. length
    
    For the given length parameter (4*f_F*L_{1-2}/D_h), find the Mach number by
    doing a solve using the secant method
    
    Parameters
    ----------
    Lparam : float
        Fanno flow length term given by (4*f_F*L_{1-2}/D_h)
    gamma : float
        Ratio of specific heats (cp/cv)
    """
    return _Ma(log(Lparam), gamma, 0.1)