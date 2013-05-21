from __future__ import division

import cython
cimport cython

import numpy as np
cimport numpy as np

from libc.math cimport M_PI as pi, log, fabs

cdef double log_RHS(double e):
    return log(e/(4*(1-e**2)**2)*(16*e**2+pi**2*(1-e**2))**0.5)

cdef _epsilon(double logWr, double epsilon0):
    
    cdef double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300
    cdef int iter = 1

    while ((iter<=3 or fabs(f)>eps) and iter<100):
        if (iter==1):
            x1 = epsilon0 
            e = x1
        if (iter==2):
            x2 = epsilon0+0.01
            e = x2
            
        if (iter>2):
            e=x2;
            
        f = log_RHS(e)-logWr;
        
        if (iter==1):
            y1=f
        if (iter>1):
            y2=f
            x3=x2-y2/(y2-y1)*(x2-x1)
            change=fabs(y2/(y2-y1)*(x2-x1))
            y1=y2; x1=x2; x2=x3
        iter += 1
    return e

def calculate_epsilon_short(logWr, epsilon0 = 0.5):
    if isinstance(logWr, np.ndarray):
        epsilon = np.zeros_like(logWr)
        for i in range(len(logWr)):
            epsilon[i] = _epsilon(logWr[i], epsilon0)
    else:
        epsilon = _epsilon(logWr, epsilon0)
        
    return epsilon
    
cdef double log_RHS_long(double e):
    return log((6*e*(pi**2-e**2*(pi**2-4))**0.5)/((2+e**2)*(1-e**2)))

cdef _epsilon_long(double logWr, double epsilon0):
    
    cdef double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,change=999,f=999,T=300
    cdef int iter = 1

    while ((iter<=3 or fabs(f)>eps) and iter<100):
        if (iter==1):
            x1 = epsilon0 
            e = x1
        if (iter==2):
            x2 = epsilon0+0.01
            e = x2
            
        if (iter>2):
            e=x2;
            
        f = log_RHS_long(e)-logWr;
        
        if (iter==1):
            y1=f
        if (iter>1):
            y2=f
            x3=x2-y2/(y2-y1)*(x2-x1)
            change=fabs(y2/(y2-y1)*(x2-x1))
            y1=y2; x1=x2; x2=x3
        iter += 1
    return e

def calculate_epsilon_long(logWr, epsilon0 = 0.1):
    if isinstance(logWr, np.ndarray):
        epsilon = np.zeros_like(logWr)
        for i in range(len(logWr)):
            epsilon[i] = _epsilon_long(logWr[i], epsilon0)
    else:
        epsilon = _epsilon_long(logWr, epsilon0)
        
    return epsilon