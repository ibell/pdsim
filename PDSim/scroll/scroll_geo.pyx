
from __future__ import division

import numpy as np
import cython

import matplotlib.pyplot as plt
from PDSim.scroll.plots import plotScrollSet
 
# A container for the values for the heat transfer angles
cdef class HTAnglesClass(object):
    def __init__(self):
        pass

#This is a list of all the members in geoVals
geoValsvarlist=['h','phi_i0','phi_is','phi_ie','phi_e',
                'phi_o0','ro','rb','phi_os','phi_oe','t',
                'xa_arc1','ya_arc1','ra_arc1','t1_arc1','t2_arc1',
                'xa_arc2','ya_arc2','ra_arc2','t1_arc2','t2_arc2',
                'b_line', 't1_line', 't2_line', 'm_line',
                'x0_wall','y0_wall','r_wall',
                'delta_radial', 'delta_flank',
                'phi_ie_offset','delta_suction_offset']
 
def rebuild_geoVals(d):
    geo = geoVals()
    for atr in geoValsvarlist:
        setattr(geo,atr,d[atr])
    return geo 
    
cdef class geoVals:
    
    def __init__(self):
        self.phi_ie_offset = 0.0
        
    def __reduce__(self):
        d={}
        for atr in geoValsvarlist:
            d[atr]=getattr(self,atr)
        return rebuild_geoVals,(d,)
    def __repr__(self):
        s='geoVals instance at '+str(id(self))+'\n'
        for atr in geoValsvarlist:
            s+=atr+': '+str(getattr(self,atr))+'\n'
        return s
        
cpdef double fxA(double rb, double phi, double phi0):
    return rb**3/3.0*(4.0*((phi-phi0)**2-2.0)*sin(phi)+(phi0-phi)*((phi-phi0)**2-8.0)*cos(phi))

cpdef double fyA(double rb, double phi, double phi0):
    return rb**3/3.0*((phi0-phi)*((phi-phi0)**2-8.0)*sin(phi)-4.0*((phi-phi0)**2-2.0)*cos(phi))

cpdef double theta_d(geoVals geo):  
    """ 
    Discharge angle
    
    Optional Parameters
    
    ======  =======================================
    key     value
    ======  =======================================
    geo     The class that defines the geometry of the compressor.
    ======  =======================================
    
    """
    N_c_max=floor((geo.phi_ie-geo.phi_os-pi)/(2*pi))
    return geo.phi_ie-geo.phi_os-2*pi*N_c_max-pi

cpdef nC_Max(geo):
    return int(floor((geo.phi_ie-geo.phi_os-pi)/(2.0*pi)))

cpdef int getNc(double theta, geoVals geo):
    """ 
    The number of pairs of compression chambers in existence at a given 
    crank angle 
    
    Arguments:
        theta : float
            The crank angle in radians.
        geo : geoVals instance

    Returns:
        Nc : int
            Number of pairs of compressions chambers
            
    """
     
    return int(floor((geo.phi_ie-theta-geo.phi_os-pi)/(2*pi)))
    
def polyarea(x,y):
    N=len(x)
    area = 0.0
    for i in range(N):
        j = (i+1) % N
        area = area + x[i]*y[j] - y[i]*x[j]
    return area/2.0
    
def polycentroid(xi,yi):
    # Add additional element if needed to close polygon
    if not xi[0]==xi[-1] or not yi[0]==yi[-1]:
        x=np.r_[xi,xi[-1]]
        y=np.r_[yi,yi[-1]]
    else:
        x=xi
        y=yi
    sumx=0.0
    sumy=0.0
    for i in range(len(x)-1):
        sumx=sumx+(x[i]+x[i+1])*(x[i]*y[i+1]-x[i+1]*y[i])
        sumy=sumy+(y[i]+y[i+1])*(x[i]*y[i+1]-x[i+1]*y[i])
    return sumx/(6*polyarea(x,y)),sumy/(6*polyarea(x,y))
    
cpdef tuple _coords_inv_np(np.ndarray[np.float_t] phi, geoVals geo,double theta, flag=""):
    """
    Internal function that does the calculation if phi is definitely a 1D numpy vector
    """
    phi_i0=geo.phi_i0
    phi_o0=geo.phi_o0
    phi_ie=geo.phi_ie
    rb=geo.rb
    ro=rb*(pi-phi_i0+phi_o0)
    om=phi_ie-theta+3.0*pi/2.0

    if flag=="fi":
        x = rb*np.cos(phi)+rb*(phi-phi_i0)*np.sin(phi)
        y = rb*np.sin(phi)-rb*(phi-phi_i0)*np.cos(phi)
    elif flag=="fo":
        x = rb*np.cos(phi)+rb*(phi-phi_o0)*np.sin(phi)
        y = rb*np.sin(phi)-rb*(phi-phi_o0)*np.cos(phi)
    elif flag=="oi":
        x = -rb*np.cos(phi)-rb*(phi-phi_i0)*np.sin(phi)+ro*np.cos(om)
        y = -rb*np.sin(phi)+rb*(phi-phi_i0)*np.cos(phi)+ro*np.sin(om)
    elif flag=="oo":
        x = -rb*np.cos(phi)-rb*(phi-phi_o0)*np.sin(phi)+ro*np.cos(om)
        y = -rb*np.sin(phi)+rb*(phi-phi_o0)*np.cos(phi)+ro*np.sin(om)
    else:
        raise ValueError('flag not valid')
    return (x,y)
    
cpdef tuple _coords_inv_d(double phi, geoVals geo,double theta, flag=""):
    """
    Internal function that does the calculation if phi is a double variable 
    """
    phi_i0=geo.phi_i0
    phi_o0=geo.phi_o0
    phi_ie=geo.phi_ie
    rb=geo.rb
    ro=rb*(pi-phi_i0+phi_o0)
    om=phi_ie-theta+3.0*pi/2.0

    if flag=="fi":
        x = rb*cos(phi)+rb*(phi-phi_i0)*sin(phi)
        y = rb*sin(phi)-rb*(phi-phi_i0)*cos(phi)
    elif flag=="fo":
        x = rb*cos(phi)+rb*(phi-phi_o0)*sin(phi)
        y = rb*sin(phi)-rb*(phi-phi_o0)*cos(phi)
    elif flag=="oi":
        x = -rb*cos(phi)-rb*(phi-phi_i0)*sin(phi)+ro*cos(om)
        y = -rb*sin(phi)+rb*(phi-phi_i0)*cos(phi)+ro*sin(om)
    elif flag=="oo":
        x = -rb*cos(phi)-rb*(phi-phi_o0)*sin(phi)+ro*cos(om)
        y = -rb*sin(phi)+rb*(phi-phi_o0)*cos(phi)+ro*sin(om)
    else:
        raise ValueError('flag not valid')
    return (x,y)    
  
cpdef tuple coords_inv(phi,geoVals geo,double theta,flag="fi"):
    """ 
    
    def coords_inv(phi,geo,theta,flag="fi")
    
    The involute angles corresponding to the points along the involutes
    (fixed inner [fi], fixed scroll outer involute [fo], orbiting
    scroll outer involute [oo], and orbiting scroll inner involute [oi] )
    
    Arguments:
        phi_vec : 1D numpy array or double
            vector of involute angles
        geo : geoVals class
            scroll compressor geometry
        theta : float
            crank angle in the range 0 to :math: `2\pi`
        flag : string
            involute of interest, possible values are 'fi','fo','oi','oo'
            
    Returns:
        (x,y) : tuple of coordinates on the scroll
    """
    if type(phi) is np.ndarray:
        return _coords_inv_np(phi,geo,theta,flag)
    else:
        return _coords_inv_d(phi,geo,theta,flag)

def coords_norm(phi_vec,geo,theta,flag="fi"):
    """ 
    The x and y coordinates of a unit normal vector pointing towards
    the scroll involute for the the involutes
    (fixed inner [fi], fixed scroll outer involute [fo], orbiting
    scroll outer involute [oo], and orbiting scroll inner involute [oi])
    
    Arguments:
        phi_vec : 1D numpy array
            vector of involute angles
        geo : geoVals class
            scroll compressor geometry
        theta : float
            crank angle in the range 0 to :math: `2\pi`
        flag : string
            involute of interest, possible values are 'fi','fo','oi','oo'
            
    Returns:
        (nx,ny) : tuple of unit normal coordinates pointing towards scroll wrap
    """
    
    phi_i0=geo.phi_i0   
    phi_o0=geo.phi_o0
    phi_ie=geo.phi_ie
    rb=geo.rb
    if not type(phi_vec) is np.ndarray:
        if type(phi_vec) is list:
            phi_vec=np.array(phi_vec)
        else:
            phi_vec=np.array([phi_vec])
    nx=np.zeros(np.size(phi_vec))
    ny=np.zeros(np.size(phi_vec))

    for i in xrange(np.size(phi_vec)):
        phi=phi_vec[i]
        if flag=="fi":
            nx[i] = +sin(phi)
            ny[i] = -cos(phi)
        elif flag=="fo":
            nx[i] = -sin(phi)
            ny[i] = +cos(phi)
        elif flag=="oi":
            nx[i] = -sin(phi)
            ny[i] = +cos(phi)
        elif flag=="oo":
            nx[i] = +sin(phi)
            ny[i] = -cos(phi)
        else:
            print "Uh oh... error in coords_norm"
    return (nx,ny)

def setDiscGeo(geo,Type='Sanden',r2=0.001,**kwargs):
    """
    Sets the discharge geometry for the compressor based on the arguments.
    Also sets the radius of the wall that contains the scroll set
    
    Arguments:
        geo : geoVals class
            class containing the scroll compressor geometry
        Type : string
            Type of discharge geometry, options are ['Sanden'],'2Arc','ArcLineArc'
        r2 : float or string
            Either the radius of the smaller arc as a float or 'PMP' for perfect meshing
            If Type is 'Sanden', this value is ignored
    
    Keyword Arguments:
    
    ========     ======================================================================
    Value        Description
    ========     ======================================================================
    r1           the radius of the large arc for the arc-line-arc solution type
    ========     ======================================================================
    
    """
    
    #Recalculate the orbiting radius
    geo.ro=geo.rb*(pi-geo.phi_i0+geo.phi_o0)
    if Type == 'Sanden':
        geo.x0_wall=0.0
        geo.y0_wall=0.0
        geo.r_wall=0.065
        setDiscGeo(geo,Type='ArcLineArc',r2=0.003178893902,r1=0.008796248080)
    elif Type == '2Arc':
        (x_is,y_is) = coords_inv(geo.phi_is,geo,0,'fi')
        (x_os,y_os) = coords_inv(geo.phi_os,geo,0,'fo')
        (nx_is,ny_is) = coords_norm(geo.phi_is,geo,0,'fi')
        (nx_os,ny_os) = coords_norm(geo.phi_os,geo,0,'fo')
        dx=x_is-x_os
        dy=y_is-y_os
        
        r2max=0
        a=cos(geo.phi_os-geo.phi_is)+1.0
        b=geo.ro*a-dx*(sin(geo.phi_os)-sin(geo.phi_is))+dy*(cos(geo.phi_os)-cos(geo.phi_is))
        c=1.0/2.0*(2.0*dx*sin(geo.phi_is)*geo.ro-2.0*dy*cos(geo.phi_is)*geo.ro-dy**2-dx**2)
        if geo.phi_os-(geo.phi_is-pi)>1e-15:
            r2max=(-b+sqrt(b**2-4.0*a*c))/(2.0*a)
        elif geo.phi_os==geo.phi_is-pi:
            r2max=-c/b
        else:
            raise ValueError('Error, must enforce phi_os > phi_is-pi to avoid scroll crashing :: phi_os %.16f phi_os %.16f, phi_is-pi %.16f' %(geo.phi_os,geo.phi_is-pi))
            
        if type(r2) is not float and r2=='PMP':
            r2=r2max
            
        if r2>r2max:
            print 'r2 is too large, max value is : %0.5f' %(r2max)
        
        xarc2 =  x_os+nx_os*r2
        yarc2 =  y_os+ny_os*r2
        
        r1=((1.0/2*dy**2+1.0/2*dx**2+r2*dx*sin(geo.phi_os)-r2*dy*cos(geo.phi_os))
               /(r2*cos(geo.phi_os-geo.phi_is)+dx*sin(geo.phi_is)-dy*cos(geo.phi_is)+r2))
        
        
        ## Negative sign since you want the outward pointing unit normal vector
        xarc1 =  x_is-nx_is*r1
        yarc1 =  y_is-ny_is*r1
                
        geo.xa_arc2=xarc2
        geo.ya_arc2=yarc2
        geo.ra_arc2=r2
        geo.t1_arc2=atan2(yarc1-yarc2,xarc1-xarc2)
        geo.t2_arc2=atan2(y_os-yarc2,x_os-xarc2)
        while geo.t2_arc2<geo.t1_arc2:
            geo.t2_arc2=geo.t2_arc2+2.0*pi;
    
        geo.xa_arc1=xarc1
        geo.ya_arc1=yarc1
        geo.ra_arc1=r1
        geo.t2_arc1=atan2(y_is-yarc1,x_is-xarc1)
        geo.t1_arc1=atan2(yarc2-yarc1,xarc2-xarc1)
        while geo.t2_arc1<geo.t1_arc1:
            geo.t2_arc1=geo.t2_arc1+2.0*pi;
        
        """ 
        line given by y=m*t+b with one element at the intersection
        point
        
        with b=0, m=y/t
        """
        geo.b_line=0.0
        geo.t1_line=xarc2+r2*cos(geo.t1_arc2)
        geo.t2_line=geo.t1_line
        geo.m_line=(yarc2+r2*sin(geo.t1_arc2))/geo.t1_line
        
        """ 
        Fit the wall to the chamber
        """
        geo.x0_wall=geo.ro/2.0*cos(geo.phi_ie-pi/2-pi)
        geo.y0_wall=geo.ro/2.0*sin(geo.phi_ie-pi/2-pi)
        (x,y)=coords_inv(geo.phi_ie,geo,pi,'fo')
        geo.r_wall=1.03*sqrt((geo.x0_wall-x)**2+(geo.y0_wall-y)**2)
    elif Type=='ArcLineArc':
        (x_is,y_is) = coords_inv(geo.phi_is,geo,0,'fi')
        (x_os,y_os) = coords_inv(geo.phi_os,geo,0,'fo')
        (nx_is,ny_is) = coords_norm(geo.phi_is,geo,0,'fi')
        (nx_os,ny_os) = coords_norm(geo.phi_os,geo,0,'fo')
        dx=x_is-x_os
        dy=y_is-y_os
        
        r2max=0
        a=cos(geo.phi_os-geo.phi_is)+1.0
        b=geo.ro*a-dx*(sin(geo.phi_os)-sin(geo.phi_is))+dy*(cos(geo.phi_os)-cos(geo.phi_is))
        c=1.0/2.0*(2.0*dx*sin(geo.phi_is)*geo.ro-2.0*dy*cos(geo.phi_is)*geo.ro-dy**2-dx**2)
        if geo.phi_os-(geo.phi_is-pi)>1e-12:
            r2max=(-b+sqrt(b**2-4.0*a*c))/(2.0*a)
        elif geo.phi_os-(geo.phi_is-pi)<1e-12:
            r2max=-c/b
        else:
            print 'error with starting angles phi_os %.16f phi_is-pi %.16f' %(geo.phi_os,geo.phi_is-pi)
            
        if type(r2) is not float and r2=='PMP':
            r2=r2max
                
        if r2>r2max:
            print 'r2 is too large, max value is : %0.5f' %(r2max)
        
        xarc2 =  x_os+nx_os*r2
        yarc2 =  y_os+ny_os*r2
        
        if 'r1' not in kwargs:
            r1=r2+geo.ro
        else:
            r1=kwargs['r1']
        
        ## Negative sign since you want the outward pointing unit normal vector
        xarc1 =  x_is-nx_is*r1
        yarc1 =  y_is-ny_is*r1
                
        geo.xa_arc2=xarc2
        geo.ya_arc2=yarc2
        geo.ra_arc2=r2
        geo.t2_arc2=atan2(y_os-yarc2,x_os-xarc2)
    
        geo.xa_arc1=xarc1
        geo.ya_arc1=yarc1
        geo.ra_arc1=r1
        geo.t2_arc1=atan2(y_is-yarc1,x_is-xarc1)
                
        alpha=atan2(yarc2-yarc1,xarc2-xarc1)
        d=sqrt((yarc2-yarc1)**2+(xarc2-xarc1)**2)
        beta=acos((r1+r2)/d)
        L=sqrt(d**2-(r1+r2)**2)
        t1=alpha+beta
        
        (xint,yint)=(xarc1+r1*cos(t1)+L*sin(t1),yarc1+r1*sin(t1)-L*cos(t1))
        t2=atan2(yint-yarc2,xint-xarc2)
        
        geo.t1_arc1=t1
#        (geo.t1_arc1,geo.t2_arc1)=sortAnglesCW(geo.t1_arc1,geo.t2_arc1)
        
        geo.t1_arc2=t2
#        (geo.t1_arc2,geo.t2_arc2)=sortAnglesCCW(geo.t1_arc2,geo.t2_arc2)

        while geo.t2_arc2<geo.t1_arc2:
            geo.t2_arc2=geo.t2_arc2+2.0*pi;
        while geo.t2_arc1<geo.t1_arc1:
            geo.t2_arc1=geo.t2_arc1+2.0*pi;
        """ 
        line given by y=m*t+b with one element at the intersection
        point
        
        with b=0, m=y/t
        """
        geo.m_line=-1/tan(t1)
        geo.t1_line=xarc1+r1*cos(geo.t1_arc1)
        geo.t2_line=xarc2+r2*cos(geo.t1_arc2)
        geo.b_line=yarc1+r1*sin(t1)-geo.m_line*geo.t1_line
        
        """ 
        Fit the wall to the chamber
        """
        geo.x0_wall=geo.ro/2.0*cos(geo.phi_ie-pi/2-pi)
        geo.y0_wall=geo.ro/2.0*sin(geo.phi_ie-pi/2-pi)
        (x,y)=coords_inv(geo.phi_ie,geo,pi,'fo')
        geo.r_wall=1.03*sqrt((geo.x0_wall-x)**2+(geo.y0_wall-y)**2)
        
    else:
        raise AttributeError('Type not understood, should be one of 2Arc or ArcLineArc')

def overlay_injection_port(theta, geo, phi, ax, inner_outer):
    """
    Plot the injection ports on an axis - no scroll wrap plot is generated.  Also see
    plot_injection_ports()
    
    Parameters
    ---------- 
    theta : float
        crank angle in the range [0, :math:`2\pi`] 
    geo : geoVals instance
    phi : float
        Involute angle in radians
    ax : matplotlib axis instance
    inner_outer : string
        If ``'i'``, phi is along the inner involute of the fixed scroll
        
        If ``'o'``, phi is along the outer involute of the fixed scroll
        
    Notes
    -----
    If you want symmetric injection ports, the ones on the inner involute 
    should have a value of phi that is pi radians greater than those on the 
    outer involute
    
    """
    
    #Common terms
    rport = geo.t/2.0
    t = np.linspace(0,2*pi,100)
    
    if inner_outer == 'o':
        #Involute angle along the outer involute of the scroll wrap
        x, y = coords_inv(phi, geo, theta, 'fo')
        nx, ny = coords_norm(phi, geo, theta, 'fo')
        xc,yc = x-nx*rport,y-ny*rport
        ax.plot(xc + rport*np.cos(t),yc+rport*np.sin(t),'k')
    elif inner_outer == 'i':
        x, y = coords_inv(phi, geo, theta, 'fi')
        nx, ny = coords_norm(phi, geo, theta, 'fi')
        xc,yc = x-nx*rport,y-ny*rport
        ax.plot(xc + rport*np.cos(t),yc+rport*np.sin(t),'k')
    else:
        raise KeyError
    
def plot_injection_ports(theta, geo, phi, ax, inner_outer):
    """
    Plot the injection ports
    
    Parameters
    ---------- 
    theta : float
        crank angle in the range [0, :math:`2\pi`] 
    geo : geoVals instance
    phi : float
        Involute angle in radians
    ax : matplotlib axis instance
    inner_outer : string
        If ``'i'``, phi is along the inner involute of the fixed scroll
        
        If ``'o'``, phi is along the outer involute of the fixed scroll
        
    Notes
    -----
    If you want symmetric injection ports, the ones on the inner involute 
    should have a value of phi that is pi radians greater than those on the 
    outer involute
    
    """
    #Plot the scrolls (symmetric)
    plotScrollSet(theta, geo, axis = ax)
    
    #Plot the port
    overlay_injection_port(theta, geo, phi, ax, inner_outer)

cpdef double min2(double a, double b):
    return a if a<b else b

cpdef double max2(double a, double b):
    return a if a>b else b
        
cpdef double radial_leakage_area(double theta, geoVals geo, bytes key1, bytes key2, location = str('up')) except *:
    """
    Get the flow area of the flow path for a given radial flow pair
    
    Parameters
    ----------
    theta : float
        crank angle in the range [0, :math:`2\pi`] 
    geo : geoVals instance
    key1 : string
    key2 : string
    location : string, one of ['up','mid','down'], optional
        What part of the wrap is used to determine the area.
    
    Returns
    -------
    Area in [\ :math:`m^2`\ ]
    
    """
    cython.declare(phi_min = cython.double, 
                   phi_max = cython.double)
    #Get the bounding angles
    phi_min,phi_max = radial_leakage_angles(theta,geo,key1,key2)
    if location =='up':
        phi_0 = geo.phi_i0
    elif location == 'down':
        phi_0 = geo.phi_o0
    elif location == 'mid':
        phi_0 = (geo.phi_i0+geo.phi_o0)/2
    else:
        raise KeyError
    A = geo.delta_radial*geo.rb*((phi_max**2-phi_min**2)/2-phi_0*(phi_max-phi_min))
    return A
    
cpdef tuple radial_leakage_angles(double theta, geoVals geo, bytes key1, bytes key2):
    """
    Get the angles for a given radial flow pair
    
    Parameters
    ----------
    theta : float
        crank angle in the range [0, :math:`2\pi`] 
    geo : geoVals instance
    key1 : string
    key2 : string
    location : string, one of ['up','mid','down'], optional
        What part of the wrap is used to determine the area.
    
    Returns
    -------
    tuple of values with (``phi_min``, ``phi_max``)
    
    Raises
    ------
    ``KeyError`` if keys are invalid or undefined at the given crank angle 
    
    """
    cython.declare(phi_min = cython.double, 
                   phi_max = cython.double, 
                   sort = cython.tuple,
                   alpha = cython.long,
                   Nc = cython.long)
    phi_min = 9e99
    phi_max = 9e99
    Nc = getNc(theta,geo)
    
    sort = tuple(sorted((key1,key2)))
    
    #These are always in existence
    if sort == ('s2','sa') or sort == ('s1','sa'):
            phi_max = geo.phi_ie
            phi_min = max2(geo.phi_ie - theta, phi_s_sa(theta,geo)+geo.phi_o0-geo.phi_i0)
    #suction chambers only in contact with each other beyond theta = pi
    elif sort == ('s1','s2'):
        if theta > pi:
            phi_max = phi_s_sa(theta,geo)+geo.phi_o0-geo.phi_i0
            phi_min = geo.phi_ie - theta
            #Ensure that at the very least phi_max is greater than  phi_min
            phi_max = max2(phi_min, phi_max)
        else:
            #They are the same so there is no flow area
            phi_max = geo.phi_ie - theta + 0.0000001
            phi_min = geo.phi_ie - theta
    
    elif Nc == 0 and phi_max>1e90:
        if sort == ('d2','s1') or sort == ('d1','s2'):
                phi_max = geo.phi_ie - theta
                phi_min = geo.phi_ie - theta - pi
        elif sort == ('d1','d2'):
                phi_max = geo.phi_ie - theta - pi
                phi_min = geo.phi_is
        elif theta > theta_d(geo):
            print 'theta = ', theta,theta_d(geo)
            print 'Nc: {Nc:d}'.format(Nc=Nc)
            raise KeyError('Nc: {Nc:d} sort {sort:s}'.format(Nc=Nc, sort = str(sort)))
    
    if Nc >= 1 and phi_max > 1e90:
        if sort == ('c2.1','sa') or sort == ('c1.1','sa'):
                phi_max = max2(geo.phi_ie - theta, phi_s_sa(theta,geo)+geo.phi_o0-geo.phi_i0 )
                phi_min = min2(geo.phi_ie - theta, phi_s_sa(theta,geo)+geo.phi_o0-geo.phi_i0 )
        elif sort == ('c2.1','s1') or sort == ('c1.1','s2'):
                #TODO: this could be improved to take into account the non-perfect separation between s-sa and phi_ie
                phi_max = geo.phi_ie - theta #this is where the change needs to be made
                phi_min = geo.phi_ie - theta - pi
        elif sort == ('c1.1','c2.1'):
                phi_max = geo.phi_ie - theta - pi
                phi_min = geo.phi_ie - theta - 2*pi
        elif Nc == 1 and (sort == ('c2.1','d1') or sort == ('c1.1','d2')):
                phi_max = geo.phi_ie - theta - 2*pi
                phi_min = geo.phi_is
        elif Nc == 1 and theta > theta_d(geo):
            print 'Nc: {Nc:d} key1: {k1:s} key2: {k2:s} theta: {theta:f} theta_d: {theta_d:f}'.format(Nc=Nc,k1=key1,k2=key2,theta = theta,theta_d = theta_d(geo))
            raise KeyError
                
    #Nc > 1
    if Nc > 1 and phi_max > 1e90: 
        for alpha in range(2, Nc+1):
            if (sort == sorted(('c2.'+str(alpha),'c1.'+str(alpha-1))) or
                sort == sorted(('c1.'+str(alpha),'c2.'+str(alpha-1)))):
                phi_max = geo.phi_ie - theta - 2*pi*(alpha-1)
                phi_min = geo.phi_ie - theta - 2*pi*(alpha-1) - pi
                break
            elif sort == sorted(('c2.'+str(alpha),'c1.'+str(alpha))):
                phi_max = geo.phi_ie - theta - 2*pi*(alpha-1) - pi
                phi_min = geo.phi_ie - theta - 2*pi*(alpha)
                break
        if phi_max > 1e90:
            if sort == ('c2.'+str(Nc),'d1') or sort == ('c1.'+str(Nc),'d2'):
                phi_max = geo.phi_ie - theta - 2*pi*Nc
                phi_min = geo.phi_is
    
    if phi_max is None or phi_min is None:
        raise KeyError ('For the pair ('+key1+','+key2+') there were no angles found')
    if phi_min > phi_max:
        raise ValueError ('For the keys ('+key1+','+key2+') @theta = '+str(theta)+' max < min (error because '+str(phi_max)+' < '+str(phi_min)+')')
    return (phi_min, phi_max)

def radial_leakage_pairs(geo):
    """
    Returns a list of all possible pairings for the radial leakages 
    
    Parameters
    ----------
    None
    
    Returns
    -------
    A list of tuples with the entries of (``key1``,``key2``) where ``key1`` and 
    ``key2`` are the keys for the control volumes
        
    Notes
    -----
    See page 125 of Bell, Ian, "Theoretical and Experimental Analysis of Liquid
    Flooded Compression in Scroll Compressors", PhD. Thesis, Purdue University,
    http://docs.lib.purdue.edu/herrick/2/
    
    Different analysis is used for 0, 1, >1 sets of compression chambers.  
    But you know from the geometry the maximum number of pairs of compression
    chambers, and can therefore determine the possible pairings *a priori*.
    
    """
    def remove_duplicates(pairs):
        #Sort each element of the list
        pairs = [sorted(pair) for pair in pairs]
        
        seen = set()
        # Keep any elements that are unique
        #
        # Sets cannot have duplicates, so adding a value to a set that is already
        # there returns a False value
        return [ x for x in pairs if str( x ) not in seen and not seen.add( str( x ) )]
    
    # You are guaranteed to have a rotational angle between 0 and 2*pi radians
    # where you lose a compression chamber            
    Nc_max = nC_Max(geo)
    Nc_min = Nc_max - 1
    
    #These are always there
    pairs = [('s1','sa'),
             ('s2','sa'),
             ('s1','s2')
             ]
    for Nc in [Nc_max, Nc_min]:
        if Nc == 0:
            pairs += [('d2','s1'),
                      ('d1','s2'),
                      ('d2','d1')]
        elif Nc == 1:
            pairs += [('c2.1','sa'),
                      ('c1.1','sa'),
                      ('c2.1','s1'),
                      ('c1.1','s2'),
                      ('c1.1','c2.1')]
            if Nc == Nc_max:
                pairs += [('d2','c1.1'),
                          ('d1','c2.1'),
                          ]
        elif Nc > 1:
            pairs += [('c2.1','sa'),
                      ('c1.1','sa'),
                      ('c2.1','s1'),
                      ('c1.1','s2'),
                      ('c2.1','c1.1'),
                      ]
            #Nc is > 1, so alpha is in the range 2, Nc inclusive
            for alpha in range(2,Nc+1):
                pairs += [('c2.'+str(alpha),'c1.'+str(alpha-1)),
                          ('c1.'+str(alpha),'c2.'+str(alpha-1)),
                          ('c1.'+str(alpha),'c2.'+str(alpha)),
                          ]
            pairs += [
                      ('d2','c1.'+str(Nc)),
                      ('d1','c2.'+str(Nc)),
                      ]
    
    return remove_duplicates(pairs)
    

cpdef HTAnglesClass HT_angles(double theta, geoVals geo, bytes key):
    """
    Return the heat transfer bounding angles for the given control volume
    
    Parameters
    ----------
    theta : float
        Crank angle in the range [:math:`0,2\pi`]
    geo : geoVals instance
    key : string
        Key for the control volume following the scroll compressor 
        naming conventions
    
    Returns
    -------
    angles : HTAngles Class 
        with the attributes:
        phi_1_i: maximum involute angle on the inner involute of the wrap 
        that forms the outer wall of the CV
        
        phi_2_i: minimum involute angle on the inner involute of the wrap 
        that forms the outer wall of the CV
        
        phi_1_o: maximum involute angle on the outer involute of the wrap 
        that forms the inner wall of the CV
        
        phi_2_o: minimum involute angle on the outer involute of the wrap 
        that forms the inner wall of the CV
        
    Notes
    -----
    The keys s1, c1.x, and d1 have as their outer wrap the fixed scroll
    
    The keys s2, c2.x, and d2 have as their outer wrap the orbiting scroll
    
    "Minimum", and "Maximum" refer to absolute values of the angles
    
    Raises
    ------
    If key is not valid, raises a KeyError
    """
    cython.declare(alpha = cython.int)
    angles = HTAnglesClass()
    ## TODO: Offset considerations to the angles
    if key == 's1' or key == 's2':
        angles.phi_1_i = geo.phi_ie
        angles.phi_2_i = geo.phi_ie-theta
        angles.phi_1_o = phi_s_sa(theta, geo)
        angles.phi_2_o = geo.phi_oe - geo.phi_o0 - pi - theta
        return angles
    elif key.startswith('c1') or key.startswith('c2'):
        alpha = int(key.split('.')[1])
        if alpha > getNc(theta,geo):
            raise KeyError('CV '+key+' does not exist')
        else:
            angles.phi_1_i = geo.phi_ie - theta - (alpha-1)*2*pi 
            angles.phi_2_i = geo.phi_ie - theta - alpha*2*pi
            angles.phi_1_o = geo.phi_oe - pi - theta - (alpha-1)*2*pi
            angles.phi_2_o = geo.phi_oe - geo.phi_o0 - pi - theta - alpha*2*pi
            return angles
    elif key == 'd1' or key == 'd2':
        alpha = getNc(theta,geo)+1
        angles.phi_1_i = geo.phi_ie - theta - geo.phi_i0 - (alpha-1)*2*pi
        angles.phi_2_i = geo.phi_is
        angles.phi_1_o = geo.phi_oe - pi - theta - (alpha-1)*2*pi
        angles.phi_2_o = geo.phi_os
        return angles
    else:
        return None

def plot_HT_angles(theta, geo, keys, involute):
    """
    Plot an involute bound angle graph for each CV for checking of  
    heat transfer bounding angles purposes.
    
    Parameters
    ----------
    theta : float
    geo : geoVals instance
    keys : list of strings of compliant CV keys
    involute : string ['i','o']
        'i': inner involute of the wrap forming the outer surface of CV
        
        'o': outer involute of the wrap forming the inner surface of CV
    
    """
    fig = plt.figure()
    ax = fig.add_axes((0.15,0.15,0.8,0.8))
    if involute == 'i':
        ax.set_title('Inner involute angles on outer surface of CV')
    elif involute == 'o':
        ax.set_title('Outer involute angles on inner surface of CV')
        
    for i, key in enumerate(keys):
        y = np.r_[i+1, i+1]
        try:
            angles = HT_angles(theta, geo, key)
            if involute == 'i':
                x = np.r_[angles.phi_2_i, angles.phi_1_i]
            elif involute == 'o':
                x = np.r_[angles.phi_2_o, angles.phi_1_o]
            ax.plot(x,y)
            ax.text(np.mean(x),y[0]+0.01,key,ha='center',va='bottom')
        except KeyError:
            pass
    ax.set_ylim(0,len(keys)+1)
    plt.show()
      
cpdef tuple SA(double theta, geoVals geo, bint poly=False, bint use_offset = True):
    """
    Volume and derivative of volume of SA chamber
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0, 2\pi`]
    geo : geoVals instance
        The geometry class
    poly : boolean
        If ``True``, also output the polygon calculation at the end of the tuple (SLOW!!)
    use_offset : boolean
        If ``True``, use the offset value from geo.phi_ie_offset, else, don't 
        use any offset
        
    Notes
    -----
    If the parameter ``geo.phi_ie_offset`` is greater than zero, the fixed scroll 
    inner involute is extended by this amount.  In general if an offset is used it
    it should be exactly :math:`\pi` radians which will bring the two suction 
    chambers together and allows for a single inlet port to the scroll set. 
    
    Returns
    -------
    V : float
    dVdTheta : float
    V_poly : float (only if ``poly = True``)
    """
    h=geo.h
    rb=geo.rb 
    phi_ie=geo.phi_ie 
    phi_o0=geo.phi_o0
    phi_oe=geo.phi_oe
    phi_i0=geo.phi_i0
    oe=geo.phi_oe
    ro=rb*(pi-phi_i0+phi_o0)
    t= rb*(phi_i0-phi_o0)
    pt = 2*pi*rb
    
    if not use_offset:
        phi_ie_offset = 0.0
    else:
        phi_ie_offset = geo.phi_ie_offset
    
    if abs(phi_ie_offset) < 1e-12:
        b=(-phi_o0+phi_ie-pi)
        D=ro/rb*((phi_i0-phi_ie)*sin(theta)-cos(theta)+1)/(phi_ie-phi_i0)
        B=1.0/2.0*(sqrt(b*b-4.0*D)-b)
        B_prime=-ro/rb*(sin(theta)+(phi_i0-phi_ie)*cos(theta))/((phi_ie-phi_i0)*sqrt(b*b-4*D))
        V_Isa=h*rb**2/6.0*(pow(phi_oe-phi_o0,3)-pow(phi_ie-pi+B-phi_o0,3))
        
        V=h*pi*geo.r_wall**2-2*V_Isa
        dV=h*rb**2*pow(phi_ie-pi+B-phi_o0,2)*B_prime
    
    else:
        """
        The suction area is a lot smaller in the case of an offset 
        scroll wrap because it is only the part around where the scroll
        wraps come close to each other.  For crank angles from 0 to pi, 
        consider the cross-section to be approximately just half the circle that
        closes the involutes.  The diameter is the pitch minus the scroll thickness
        
        Beyond theta = pi, the area grows as the s1 chamber begins to compress
        The area remaining in the chamber at theta = 2*pi becomes the s1 chamber 
        """ 
        
        Dcircle = pt-t
        Vcircle = 0.5 * (pi * Dcircle**2 / 4) * h
        dVcircle = 0.0
        
        if theta <= pi:
            V = Vcircle
            dV = dVcircle
        else:
            #Using the old style of definition
            b = (-phi_o0+phi_ie+phi_ie_offset-pi)    
            B = -ro/rb*sin(theta)/b
            B_prime = -ro/rb/b*cos(theta)
            
#            print 'theta',theta 
            #The scroll wrap portion of the outer wall of the chamber
            VO = h*rb**2/6.0*((phi_ie+phi_ie_offset-phi_i0)**3-(phi_ie+phi_ie_offset-(theta-pi)-phi_i0)**3)
            dVO = h*rb**2/2.0*((phi_ie+phi_ie_offset-(theta-pi)-phi_i0)**2)
            
#            print 'O phi',phi_ie+phi_ie_offset,phi_ie+phi_ie_offset-(theta-pi)
            # From the fixed origin to the fixed scroll using the same bounding angles
            # as for the orbiting scroll
            VIa=h*rb**2/6.0*((phi_ie-phi_o0)**3-(phi_ie-(theta-pi)-phi_o0)**3)
            dVIa=h*rb**2/2.0*((phi_ie-(theta-pi)-phi_o0)**2)
            
#            print 'I phi',phi_ie,phi_ie-(theta-pi)
            
            VIb=h*rb*ro/2.0*((phi_ie-pi+B+phi_ie_offset-phi_o0)*sin(B+phi_ie_offset+theta)+cos(B+phi_ie_offset+theta))
            dVIb=h*rb*ro*(B_prime+1)/2.0*((phi_ie-pi+B+phi_ie_offset-phi_o0)*cos(B+phi_ie_offset+theta)-sin(B+phi_ie_offset+theta))
                
            VIc=h*rb*ro/2.0
            dVIc=0.0
                
            VI = VIa+VIb-VIc
            dVI = dVIa+dVIb-dVIc
        
            V = Vcircle + VO - VI
            dV = dVcircle + dVO - dVI
         
    if not poly:
        return V,dV
    else:
        raise ValueError('no polygons for volumes of SA yet')
    
cpdef dict SA_forces(double theta, geoVals geo, bint poly = False, bint use_offset = False):
    
    h=geo.h
    rb=geo.rb 
    phi_ie=geo.phi_ie 
    phi_o0=geo.phi_o0
    phi_oe=geo.phi_oe
    phi_i0=geo.phi_i0
    ro=rb*(pi-phi_i0+phi_o0)
    t= rb*(phi_i0-phi_o0)

    phi_ie_offset = geo.phi_ie_offset
    
    if phi_ie_offset > 1e-12:
        # Calculations for break angle seem to be messed up, rolling back to the
        # older method which seems to work pretty well, though imperfectly.
        # Good enough for now
        
        b = (-phi_o0+phi_ie+phi_ie_offset-pi)    
        B = -ro/rb*sin(theta)/b
        B_prime = -ro/rb/b*cos(theta) 
        
    else:
        b=(-phi_o0+phi_ie-pi)
        D=ro/rb*((phi_i0-phi_ie)*sin(theta)-cos(theta)+1)/(phi_ie-phi_i0)
        B=1.0/2.0*(-b+sqrt(b**2-4.0*D))
        B_prime=-ro/rb*(sin(theta)+(phi_i0-phi_ie)*cos(theta))/((phi_ie-phi_i0)*sqrt(b**2-4*D))
        
    fx_p=rb*h*(cos(phi_ie)*phi_o0+sin(phi_ie)-phi_ie*cos(phi_ie) - sin(phi_ie-pi+B)-(phi_o0-phi_ie+pi-B)*cos(phi_ie-pi+B) )
    fy_p=-rb*h*((phi_ie-phi_o0)*sin(phi_ie)+cos(phi_ie) - cos(phi_ie-pi+B) - (phi_ie-pi-phi_o0+B)*sin(phi_ie-pi+B))
    M_O_p=-h*rb**2*(B-2*phi_o0+2*phi_ie-2*pi)*(B-pi)/2
    
    exact_dict = dict(fx_p = fx_p,
                      fy_p = fy_p,
                      fz_p = 0,
                      M_O_p = M_O_p,
                      cx = 0,
                      cy = 0
                      )
    
#    ############### Polygon calculations ##################
#    phi=np.linspace(phi_ie-pi+B,phi_ie,2000)
#    (xo,yo)=coords_inv(phi, geo, theta, 'oo')
#    ############### Numerical Force Calculations ###########
#    phi=np.linspace(phi_ie-pi+B,phi_ie,2000)
#    nx=np.zeros_like(phi)
#    ny=np.zeros_like(phi)
#    (nx,ny)=coords_norm(phi,geo,theta,'oo')
#    L=len(xo)
#    dA=h*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
#    dfxp_poly=dA*(nx[1:L]+nx[0:L-1])/2.0
#    dfyp_poly=dA*(ny[1:L]+ny[0:L-1])/2.0
#    fxp_poly=np.sum(dfxp_poly)
#    fyp_poly=np.sum(dfyp_poly)
#    
#    print 'fx_p',fx_p,fxp_poly
#    print 'fy_p',fy_p,fyp_poly
    
    return exact_dict

cpdef tuple S1(double theta, geoVals geo, bint poly=False, double theta_0_volume=1e-9, bint use_offset = True):
    """
    Volume and derivative of volume of S1 chamber
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0, 2\pi`]
    geo : geoVals instance
        The geometry class
    poly : boolean
        If ``True``, also output the polygon calculation at the end of the tuple (SLOW!!)
    use_offset : boolean
        If ``True``, use the offset value from geo.phi_ie_offset, else, don't 
        use any offset
        
    Notes
    -----
    If the parameter ``geo.phi_ie_offset`` is greater than zero, the fixed scroll 
    inner involute is extended by this amount.  In general if an offset is used it
    it should be exactly :math:`\pi` radians which will bring the two suction 
    chambers together and allows for a single inlet port to the scroll set. 
    
    Returns
    -------
    V : float
    dVdTheta : float
    V_poly : float (only if ``poly = True``)
    """
        
    h=geo.h
    rb=geo.rb 
    phi_ie=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    ro=rb*(pi-phi_i0+phi_o0)

    if not use_offset:
        phi_ie_offset = 0.0
    else:
        phi_ie_offset = geo.phi_ie_offset
        
    if phi_ie_offset > 1e-12:
        # Calculations for break angle seem to be messed up, rolling back to the
        # older method which seems to work pretty well, though imperfectly.
        # Good enough for now
        
        b = (-phi_o0+phi_ie+phi_ie_offset-pi)    
        B = -ro/rb*sin(theta)/b
        B_prime = -ro/rb/b*cos(theta) 
        
    else:
        b=(-phi_o0+phi_ie-pi)
        D=ro/rb*((phi_i0-phi_ie)*sin(theta)-cos(theta)+1)/(phi_ie-phi_i0)
        B=1.0/2.0*(-b+sqrt(b**2-4.0*D))
        B_prime=-ro/rb*(sin(theta)+(phi_i0-phi_ie)*cos(theta))/((phi_ie-phi_i0)*sqrt(b**2-4*D))
    
    if phi_ie_offset > 0 and theta >= pi:
        
        # Beyond a crank angle of pi radians, the offset chamber is "sealed off" 
        # besides any leakage through the channel, and don't need to find
        # the bounding involute angle on the opposing scroll wrap.
        # 
        # You can just use the compression chamber volume relationship - a linear
        # decrease in volume with the crank angle.  The compression chamber
        # index alpha is zero because this is not a conventional compression chamber
        # but a weird semi-compression chamber formed by the offset pocket
        
        Vs = -pi*h*rb*ro*(2*theta-2*phi_ie-pi+phi_i0+phi_o0)
        dVs = -2.0*pi*h*rb*ro
    
    else:
        
        # This block of code is used whenever 
    
        VO=h*rb**2/6.0*((phi_ie+phi_ie_offset-phi_i0)**3-(phi_ie-theta-phi_i0)**3)
        dVO=h*rb**2/2.0*((phi_ie-theta-phi_i0)**2)
        
        VIa=h*rb**2/6.0*((phi_ie+phi_ie_offset-pi+B-phi_o0)**3-(phi_ie-pi-theta-phi_o0)**3)
        dVIa=h*rb**2/2.0*((phi_ie+phi_ie_offset-pi+B-phi_o0)**2*B_prime+(phi_ie-pi-theta-phi_o0)**2)
            
        VIb=h*rb*ro/2.0*((phi_ie-pi+B+phi_ie_offset-phi_o0)*sin(B+phi_ie_offset+theta)+cos(B+phi_ie_offset+theta))
        dVIb=h*rb*ro*(B_prime+1)/2.0*((phi_ie-pi+B+phi_ie_offset-phi_o0)*cos(B+phi_ie_offset+theta)-sin(B+phi_ie_offset+theta))
            
        VIc=h*rb*ro/2.0
        dVIc=0.0
            
        VI = VIa+VIb-VIc
        dVI = dVIa+dVIb-dVIc
        
        Vs=VO-VI
        dVs=dVO-dVI
    
    # Add on the ficticious volume to correct for there actually being no volume
    # at theta=0
    Vs+=theta_0_volume
    
    if poly==False:
        return Vs,dVs
    elif poly==True:
        
        ############### Polygon calculations ##################
        phi=np.linspace(phi_ie-theta,phi_ie+phi_ie_offset,2000)
        (xi,yi)=coords_inv(phi, geo, theta, 'fi')
        phi=np.linspace(phi_ie-pi+B+phi_ie_offset,phi_ie-pi-theta,2000)
        (xo,yo)=coords_inv(phi, geo, theta, 'oo')
        V_poly=h*polyarea(np.r_[xi,xo,xi[0]], np.r_[yi,yo,yi[0]])

        #################################################################
        ######### Polygon plotting code for debugging purposes ##########
        #################################################################
#        if theta <= pi:
#            import pylab
#            plotScrollSet(theta,geo,shaveOn=False,offsetScroll=phi_ie_offset>0)
#            pylab.gca().fill(np.r_[xi,xo,xi[0]], np.r_[yi,yo,yi[0]])
#            pylab.gca().plot(0,0,'x')
#            pylab.gca().plot(ro*cos(phi_ie-pi/2-theta),ro*sin(phi_ie-pi/2-theta),'o')
#            
#            phi=np.linspace(phi_ie+phi_ie_offset,phi_ie-theta,2000)
#            (xo,yo)=coords_inv(phi, geo, theta, 'fi')
#            pylab.gca().plot(xo,yo,'-',lw=2)
#            pylab.gca().fill(np.r_[xo[::-1],0], np.r_[yo[::-1],0],'k',alpha = 0.5)
#            print 'VO',h*polyarea(np.r_[xo[::-1],0], np.r_[yo[::-1],0]),VO
#            
#            #VI calc
#            phi=np.linspace(phi_ie+phi_ie_offset-pi+B,phi_ie-pi-theta,2000)
#            (xo,yo)=coords_inv(phi, geo, theta, 'oo')
#            pylab.gca().plot(xo,yo,'-',lw=2)
#            pylab.gca().fill(np.r_[xo[::-1],0], np.r_[yo[::-1],0],'g',alpha = 0.5)
#            print 'VI', h*polyarea(np.r_[xo[::-1],0], np.r_[yo[::-1],0]), VI
#            VI_poly = h*polyarea(np.r_[xo[::-1],0], np.r_[yo[::-1],0])
#            
#            phi=np.linspace(phi_ie+phi_ie_offset-pi+B,phi_ie-pi-theta,2000)
#            (xo,yo)=coords_inv(phi, geo, theta, 'fo')
#            pylab.gca().fill(np.r_[xo[::-1],0], np.r_[yo[::-1],0],'b',alpha = 0.3,lw=5)
#            print 'VIa', h*polyarea(np.r_[xo[::-1],0], np.r_[yo[::-1],0]), VIa
#            
#            x1,y1 = 0.0,0.0
#            x2,y2 = coords_inv(phi_ie-pi+B+phi_ie_offset, geo, theta, 'fo')
#            x3,y3 = ro*cos(phi_ie-pi/2-theta),ro*sin(phi_ie-pi/2-theta)
#            pylab.gca().fill(np.r_[x1,x2,x3,x1], np.r_[y1,y2,y3,y1],'r',alpha = 0.8)
#            print 'VIb',h*polyarea(np.r_[x1,x2,x3,x1], np.r_[y1,y2,y3,y1]),VIb
#            
#            x1,y1 = 0.0,0.0
#            x2,y2 = coords_inv(phi_ie-pi-theta, geo, theta, 'fo')
#            x3,y3 = ro*cos(phi_ie-pi/2-theta),ro*sin(phi_ie-pi/2-theta)
#            pylab.gca().fill(np.r_[x1,x2,x3,x1], np.r_[y1,y2,y3,y1],'y',alpha=0.8)
#            print 'VIc',h*polyarea(np.r_[x1,x2,x3,x1], np.r_[y1,y2,y3,y1]), VIc
#            
#            pylab.show()
            
        return Vs,dVs,V_poly

cpdef dict S1_forces(double theta, geoVals geo, bint poly = False, double theta_0_volume=1e-9, bint use_offset = True):
    
    import warnings
    if geo.phi_ie_offset>1e-12:
        warnings.warn('S1_forces not fixed for offset scroll', RuntimeWarning)
        
    h=geo.h
    rb=geo.rb 
    phi_ie=geo.phi_ie
    phi_e=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    ro=rb*(pi-phi_i0+phi_o0)
    
    b=(-phi_o0+phi_e-pi)
    D=ro/rb*((phi_i0-phi_e)*sin(theta)-cos(theta)+1)/(phi_e-phi_i0)
    B=1.0/2.0*(sqrt(b**2-4.0*D)-b)
    B_prime=-ro/rb*(sin(theta)+(phi_i0-phi_ie)*cos(theta))/((phi_e-phi_i0)*sqrt(b**2-4*D))
    
    VO=h*rb**2/6.0*((phi_e-phi_i0)**3-(phi_e-theta-phi_i0)**3)
    dVO=h*rb**2/2.0*((phi_e-theta-phi_i0)**2)
    
    VIa=h*rb**2/6.0*((phi_e-pi+B-phi_o0)**3-(phi_e-pi-theta-phi_o0)**3)
    dVIa=h*rb**2/2.0*((phi_e-pi+B-phi_o0)**2*B_prime+(phi_e-pi-theta-phi_o0)**2)
        
    VIb=h*rb*ro/2.0*((B-phi_o0+phi_e-pi)*sin(B+theta)+cos(B+theta))
    dVIb=h*rb*ro*(B_prime+1)/2.0*((phi_e-pi+B-phi_o0)*cos(B+theta)-sin(B+theta))
    
    VIc=h*rb*ro/2.0
    dVIc=0.0
    
    #Add the small volume to avoid division by zero at theta=0
    VO += theta_0_volume
    VIa += theta_0_volume
    VIb += theta_0_volume
        
    Vs=VO-(VIa+VIb-VIc)
    dVs=dVO-(dVIa+dVIb-dVIc)
    
    cx_O=h/VO*(fxA(rb,phi_ie,phi_i0)-fxA(rb,phi_ie-theta,phi_i0))
    cy_O=h/VO*(fyA(rb,phi_ie,phi_i0)-fyA(rb,phi_ie-theta,phi_i0))
    
    cx_Ia=h/VIa*(fxA(rb,phi_ie-pi+B,phi_o0)-fxA(rb,phi_ie-pi-theta,phi_o0))
    cy_Ia=h/VIa*(fyA(rb,phi_ie-pi+B,phi_o0)-fyA(rb,phi_ie-pi-theta,phi_o0))
    
    cx_Ib=1.0/3.0*(-rb*(B-phi_o0+phi_e-pi)*sin(B+phi_e)-rb*cos(B+phi_e)-ro*sin(theta-phi_e))
    cy_Ib=1.0/3.0*(-rb*sin(B+phi_e)+rb*(B-phi_o0+phi_e-pi)*cos(B+phi_e)-ro*cos(theta-phi_e))
    
    cx_Ic=1.0/3.0*(rb*(-theta-phi_o0+phi_e-pi)*sin(theta-phi_e)-ro*sin(theta-phi_e)-rb*cos(theta-phi_e))
    cy_Ic=1.0/3.0*(rb*sin(theta-phi_e)+rb*(-theta-phi_o0+phi_e-pi)*cos(theta-phi_e)-ro*cos(theta-phi_e))

    cx_I=-(cx_Ia*VIa+cx_Ib*VIb-cx_Ic*VIc)/(VIa+VIb-VIc)+ro*cos(phi_ie-pi/2.0-theta)
    cy_I=-(cy_Ia*VIa+cy_Ib*VIb-cy_Ic*VIc)/(VIa+VIb-VIc)+ro*sin(phi_ie-pi/2.0-theta)
    
    cx=(cx_O*VO-cx_I*(VIa+VIb-VIc))/Vs
    cy=(cy_O*VO-cy_I*(VIa+VIb-VIc))/Vs
    
    if not use_offset:
        fx_p = -rb*h*(sin(B+phi_e)-(B-phi_o0+phi_e-pi)*cos(B+phi_e)+sin(theta-phi_e)-(theta+phi_o0-phi_e+pi)*cos(theta-phi_e))
        fy_p = rb*h*((B-phi_o0+phi_e-pi)*sin(B+phi_e)+cos(B+phi_e)-(theta+phi_o0-phi_e+pi)*sin(theta-phi_e)-cos(theta-phi_e))
    else:
        fx_p = -rb*h*(sin(B+phi_e)-(B-phi_o0+phi_e-pi)*cos(B+phi_e)+sin(theta-phi_e)-(theta+phi_o0-phi_e+pi)*cos(theta-phi_e))
        fy_p = rb*h*((B-phi_o0+phi_e-pi)*sin(B+phi_e)+cos(B+phi_e)-(theta+phi_o0-phi_e+pi)*sin(theta-phi_e)-cos(theta-phi_e))
        
    M_O_p=(h*rb**2*(B-theta-2*phi_o0+2*phi_e-2*pi)*(B+theta))/2
    fz_p = Vs/h
    exact_dict = dict(fx_p = fx_p,
                      fy_p = fy_p,
                      fz_p = fz_p,
                      M_O_p = M_O_p,
                      cx = cx,
                      cy = cy
                      )
    if not poly:
        return exact_dict
    else:
        ############### Polygon calculations ##################
        phi=np.linspace(phi_ie-theta,phi_ie,2000)
        (xi,yi)=coords_inv(phi, geo, theta, 'fi')
        phi=np.linspace(phi_ie-pi+B,phi_ie-pi-theta,2000)
        (xo,yo)=coords_inv(phi, geo, theta, 'oo')
        V_poly=h*polyarea(np.r_[xi,xo,xi[0]], np.r_[yi,yo,yi[0]])
        if V_poly>0.0:
            (cx_poly,cy_poly)=polycentroid(np.r_[xi,xo,xi[0]], np.r_[yi,yo,yi[0]])
        else:
            (cx_poly,cy_poly)=(xi[0],yi[0])
        ############### Numerical Force Calculations ###########
        phi=np.linspace(phi_ie-pi+B,phi_ie-pi-theta,2000)
        nx=np.zeros_like(phi)
        ny=np.zeros_like(phi)
        (nx,ny)=coords_norm(phi,geo,theta,'oo')
        L=len(xo)
        dA=h*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
        dfxp_poly=dA*(nx[1:L]+nx[0:L-1])/2.0
        dfyp_poly=dA*(ny[1:L]+ny[0:L-1])/2.0
        fxp_poly=np.sum(dfxp_poly)
        fyp_poly=np.sum(dfyp_poly)
        rOx=xo-geo.ro*cos(phi_e-pi/2-theta)
        rOx=(rOx[1:L]+rOx[0:L-1])/2
        rOy=yo-geo.ro*sin(phi_e-pi/2-theta)
        rOy=(rOy[1:L]+rOy[0:L-1])/2
        MO_poly=np.sum(rOx*dfyp_poly-rOy*dfxp_poly)
        poly_dict = dict(MO_poly = MO_poly,
                         fxp_poly = fxp_poly,
                         fyp_poly = fyp_poly,
                         cx_poly = cx_poly,
                         cy_poly = cy_poly
                         )
        exact_dict.update(poly_dict)
        return exact_dict
    
cpdef tuple S2(double theta, geoVals geo, bint poly = False, double theta_0_volume = 1e-9):
    """
    Volume and derivative of volume of S2 chamber
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0, 2\pi`]
    geo : geoVals instance
        The geometry class
    poly : boolean
        If true, also output the polygon calculation at the end of the tuple (SLOW!!)
        
    Returns
    -------
    V : float
    dVdTheta : float
    V_poly : float (only if ``poly = True``)
    """
    return S1(theta, geo, theta_0_volume = theta_0_volume, use_offset = False)
        
cpdef dict S2_forces(double theta, geoVals geo, bint poly=False, double theta_0_volume = 1e-9):
    """
    Force terms for S2 chamber
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0,2\pi`]
    geo : geoVals instance
        The geometry class
    poly : boolean, optional
        If true, also output the polygon calculations to the dict (SLOW!!)
    
    Returns
    -------
    values : dictionary
        A dictionary with fields for the analytic and numerical solutions (if requested)
    
    """
    cython.declare(h = cython.double, 
                   rb = cython.double, 
                   phi_ie = cython.double,
                   phi_o0 = cython.double,
                   phi_i0 = cython.double,
                   cx_s1 = cython.double,
                   cy_s1 = cython.double,
                   cx = cython.double,
                   cy = cython.double,
                   exact_dict = cython.dict
                   )
    h=geo.h
    rb=geo.rb
    phi_ie=geo.phi_ie
    phi_e=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    ro=rb*(pi-phi_i0+phi_o0)
    
    S1_terms = S1_forces(theta,geo, theta_0_volume = theta_0_volume, use_offset = False)
    cx_s1 = S1_terms['cx']
    cy_s1 = S1_terms['cy']
    
    (cx,cy)=(-cx_s1+ro*cos(phi_ie-pi/2-theta),-cy_s1+ro*sin(phi_ie-pi/2-theta))

    fx_p=-rb*h*(sin(theta-phi_e)-(theta+phi_i0-phi_e)*cos(theta-phi_e)+cos(phi_e)*(phi_i0-phi_e)+sin(phi_e))
    fy_p=-rb*h*((theta+phi_i0-phi_e)*sin(theta-phi_e)+cos(theta-phi_e)+sin(phi_e)*(phi_i0-phi_e)-cos(phi_e))
    M_O_p=(h*rb**2*theta*(theta+2*phi_i0-2*phi_e))/2
    fz_p = S1_terms['fz_p'] #By symmetry
    
    exact_dict = dict(fx_p = fx_p,
                      fy_p = fy_p,
                      fz_p = fz_p,
                      M_O_p = M_O_p,
                      cx = cx,
                      cy = cy
                      )
    
    if not poly:
        return exact_dict
    else:
        raise NotImplementedError('S2_forces polygon not implemented')
#        ############### Numerical Force Calculations ###########
#        phi=np.linspace(phi_ie-theta,phi_ie,2000)
#        (xo,yo)=coords_inv(phi, geo, theta, 'oi')
#        nx=np.zeros_like(phi)
#        ny=np.zeros_like(phi)
#        (nx,ny)=coords_norm(phi,geo,theta,'oi')
#        L=len(xo)
#        dA=h*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
#        fxp_poly=np.sum(dA*(nx[1:L]+nx[0:L-1])/2.0)
#        fyp_poly=np.sum(dA*(ny[1:L]+ny[0:L-1])/2.0)
#        poly_dict = dict(fxp_poly = fxp_poly,
#                         fyp_poly = fyp_poly,
##                         MO_poly = MO_poly,
##                         cx_poly = cx_poly,
##                         cy_poly = cy_poly
#                         )
#        exact_dict.update(poly_dict)
#        return exact_dict

cpdef tuple C1(double theta, int alpha, geoVals geo, bint poly=False):
    """
    Volume terms for C1,alpha chamber
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0,2\pi`]
    alpha : int
        The index of the compression chamber ( 1=outermost chamber )
    geo : geoVals instance
        The geometry class
    poly : boolean, optional
        If true, also output the polygon calculations to the dict (SLOW!!)
    
    Returns
    -------
    values : tuple
        A tuple with volume,derivative of volume and volume from polygon(if requested)
    """
    h=geo.h
    rb=geo.rb
    phi_ie=geo.phi_ie
    phi_e=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    ro=rb*(pi-phi_i0+phi_o0)
    
    ##################### Analytic Calculations ####################
    V=-pi*h*rb*ro*(2*theta+4*alpha*pi-2*phi_ie-pi+phi_i0+phi_o0)
    dV=-2.0*pi*h*rb*ro
    
    if not poly:
        return V,dV
    else:
        ##################### Polygon Calculations #####################
        phi=np.linspace(geo.phi_ie-theta-2*pi*alpha,geo.phi_ie-theta-2*pi*(alpha-1), 1000)
        (xi,yi)=coords_inv(phi, geo, theta, 'fi')
        phi=np.linspace( geo.phi_ie-theta-2*pi*(alpha-1)-pi,geo.phi_ie-theta-2*pi*alpha-pi,1000)
        (xo,yo)=coords_inv(phi, geo, theta, 'oo')
        V_poly=h*polyarea(np.r_[xi,xo], np.r_[yi,yo])
        return V,dV,V_poly

cpdef dict C1_forces(double theta, int alpha, geoVals geo, bint poly = False):
    """
    Force terms for C1,alpha chamber
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0,2\pi`]
    alpha : int
        The index of the compression chamber ( 1=outermost chamber )
    geo : geoVals instance
        The geometry class
    poly : boolean, optional
        If true, also output the polygon calculations to the dict (SLOW!!)
    
    Returns
    -------
    values : dictionary
        A dictionary with fields for the analytic and numerical solutions (if requested)
    
    """
    h=geo.h
    rb=geo.rb
    phi_ie=geo.phi_ie
    phi_e=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    ro=rb*(pi-phi_i0+phi_o0)
    
    psi=rb/3.0*(3.0*theta**2+6.0*phi_o0*theta+3.0*phi_o0**2+pi**2-15.0+(theta+phi_o0)*(12.0*pi*alpha-6.0*phi_ie)+3.0*phi_ie**2+12.0*pi*alpha*(pi*alpha-phi_ie))/(2.0*theta+phi_o0-2.0*phi_ie+phi_i0+4.0*pi*alpha-pi)
    cx=-2.0*rb*cos(theta-phi_ie)-psi*sin(theta-phi_ie)
    cy=+2.0*rb*sin(theta-phi_ie)-psi*cos(theta-phi_ie) 
    fx_p= 2.0*pi*rb*h*cos(theta-phi_e)
    fy_p=-2.0*pi*rb*h*sin(theta-phi_e)
    M_O_p=-2*pi*h*rb*rb*(theta+phi_o0-phi_e+2*pi*alpha)
    fz_p = C1(theta,alpha,geo)[0]/h
    exact_dict = dict(fx_p = fx_p,
                      fy_p = fy_p,
                      fz_p = fz_p,
                      M_O_p = M_O_p,
                      cx = cx,
                      cy = cy
                      )
    
    if not poly:
        return exact_dict
    else:
         ##################### Polygon Calculations #####################
        phi=np.linspace(geo.phi_ie-theta-2*pi*alpha,geo.phi_ie-theta-2*pi*(alpha-1), 1000)
        (xi,yi)=coords_inv(phi, geo, theta, 'fi')
        phi=np.linspace( geo.phi_ie-theta-2*pi*(alpha-1)-pi,geo.phi_ie-theta-2*pi*alpha-pi,1000)
        (xo,yo)=coords_inv(phi, geo, theta, 'oo')
        V_poly=h*polyarea(np.r_[xi,xo], np.r_[yi,yo])
        (cx_poly,cy_poly)=polycentroid(np.r_[xi,xo], np.r_[yi,yo])
        ##################### Force Calculations #########################
        phi=np.linspace( geo.phi_ie-theta-2*pi*(alpha)-pi,geo.phi_ie-theta-2*pi*(alpha-1)-pi,1000)
        (xo,yo)=coords_inv(phi, geo, theta, 'oo')
        nx=np.zeros_like(phi)
        ny=np.zeros_like(phi)
        (nx,ny)=coords_norm(phi,geo,theta,'oo')
        L=len(xo)
        dA=h*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
        dfxp_poly=dA*(nx[1:L]+nx[0:L-1])/2.0
        dfyp_poly=dA*(ny[1:L]+ny[0:L-1])/2.0
        fxp_poly=np.sum(dfxp_poly)
        fyp_poly=np.sum(dfyp_poly)
        rOx=xo-geo.ro*cos(phi_e-pi/2-theta)
        rOx=(rOx[1:L]+rOx[0:L-1])/2
        rOy=yo-geo.ro*sin(phi_e-pi/2-theta)
        rOy=(rOy[1:L]+rOy[0:L-1])/2
        MO_poly=np.sum(rOx*dfyp_poly-rOy*dfxp_poly)
        poly_dict = dict(fxp_poly = fxp_poly,
                         fyp_poly = fyp_poly,
                         MO_poly = MO_poly,
                         cx_poly = cx_poly,
                         cy_poly = cy_poly
                         )
        exact_dict.update(poly_dict)
        return exact_dict
    
cpdef tuple C2(double theta, int alpha, geoVals geo, bint poly=False):
    """
    Volume terms for C2,alpha chamber
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0,2\pi`]
    alpha : int
        The index of the compression chamber ( 1=outermost chamber )
    geo : geoVals instance
        The geometry class
    poly : boolean, optional
        If true, also output the polygon calculations to the dict (SLOW!!)
    
    Returns
    -------
    values : tuple
        A tuple with volume,derivative of volume and volume from polygon(if requested)
    """
    
    #Use the symmetry - chambers have the same volumes and derivative of volume
    return C1(theta,alpha,geo,poly)
    
cpdef dict C2_forces(double theta, int alpha, geoVals geo, bint poly=False):
    """
    Force terms for C2,alpha chamber
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0,2\pi`]
    alpha : int
        The index of the compression chamber ( 1=outermost chamber )
    geo : geoVals instance
        The geometry class
    poly : boolean, optional
        If true, also output the polygon calculations to the dict (SLOW!!)
    
    Returns
    -------
    values : dictionary
        A dictionary with fields for the analytic and numerical solutions (if requested)
    """

    h=geo.h
    rb=geo.rb
    phi_ie=geo.phi_ie
    phi_e=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    ro=geo.ro
    
    C1_dict = C1_forces(theta,alpha,geo,poly)
    cxc1 = C1_dict['cx']
    cyc1 = C1_dict['cy']
    (cx,cy)=(-cxc1+ro*cos(phi_ie-pi/2-theta),-cyc1+ro*sin(phi_ie-pi/2-theta))
    fx_p= 2.0*pi*rb*h*cos(theta-phi_e)
    fy_p=-2.0*pi*rb*h*sin(theta-phi_e)
    M_O_p=2*pi*h*rb*rb*(theta+phi_i0-phi_e+2*pi*alpha-pi)
    fz_p = C2(theta,alpha,geo)[0]/h 
    
    exact_dict = dict(fx_p = fx_p,
                      fy_p = fy_p,
                      fz_p = fz_p,
                      M_O_p = M_O_p,
                      cx = cx,
                      cy = cy
                      )
    
    if not poly:
        return exact_dict
    else:
        raise NotImplementedError('C2_forces polygon not implemented')
#        ##################### Force Calculations #########################
#        phi=np.linspace( geo.phi_ie-theta-2*pi*(alpha),geo.phi_ie-theta-2*pi*(alpha-1),1000)
#        (xo,yo)=coords_inv(phi, geo, theta, 'oi')
#        nx=np.zeros_like(phi)
#        ny=np.zeros_like(phi)
#        (nx,ny)=coords_norm(phi,geo,theta,'oi')
#        L=len(xo)
#        
#        dA=h*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
#        fxp_poly=np.sum(dA*(nx[1:L]+nx[0:L-1])/2.0)
#        fyp_poly=np.sum(dA*(ny[1:L]+ny[0:L-1])/2.0)
#        poly_dict = dict(fxp_poly = fxp_poly,
#                         fyp_poly = fyp_poly,
##                         MO_poly = MO_poly,
##                         cx_poly = cx_poly,
##                         cy_poly = cy_poly
#                         )
#        exact_dict.update(poly_dict)
#        return exact_dict

cpdef tuple D1(double theta, geoVals geo, bint poly = False):
    """
    Volume terms for D1
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0,2\pi`]
    geo : geoVals instance
        The geometry class
    poly : boolean, optional
        If true, also output the polygon calculations to the dict (SLOW!!)
    
    Returns
    -------
    values : tuple
        A tuple with volume,derivative of volume and volume from polygon(if requested)
    """
    cython.declare(Nc = cython.double)
    
    hs=geo.h
    rb=geo.rb
    phi_ie=geo.phi_ie
    phi_e=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    phi_is=geo.phi_is
    phi_os=geo.phi_os
    ro=rb*(pi-phi_i0+phi_o0)
    Nc=getNc(theta, geo)
    
    phi2=phi_ie-theta-2.0*pi*Nc
    phi1=phi_os+pi
    VO=hs*rb**2/6.0*((phi2-phi_i0)**3-(phi1-phi_i0)**3)
    dVO=-hs*rb**2/2.0*((phi2-phi_i0)**2)
    
    phi2=phi_ie-theta-2.0*pi*Nc-pi
    phi1=phi_os
    VIa=hs*rb**2/6.0*((phi2-phi_o0)**3-(phi1-phi_o0)**3)
    dVIa=-hs*rb**2/2.0*((phi2-phi_o0)**2)
    
    VIb=hs*rb*ro/2.0*((phi_os-phi_o0)*sin(theta+phi_os-phi_ie)+cos(theta+phi_os-phi_ie))
    dVIb=hs*rb*ro/2.0*((phi_os-phi_o0)*cos(theta+phi_os-phi_ie)-sin(theta+phi_os-phi_ie))
    
    VIc=hs*rb*ro/2.0
    dVIc=0.0
    
    VId= hs*rb*ro/2.0*((phi_os-phi_i0+pi)*sin(theta+phi_os-phi_ie)+cos(theta+phi_os-phi_ie)+1)
    dVId=hs*rb*ro/2.0*((phi_os-phi_i0+pi)*cos(theta+phi_os-phi_ie)-sin(theta+phi_os-phi_ie))
    
    VI=VIa+VIb+VIc+VId
    dVI=dVIa+dVIb+dVIc+dVId
    
    Vd1=VO-VI
    dVd1=dVO-dVI
    
    if not poly:
        return Vd1,dVd1
    else:
        ######################### Polygon calculations ##################
        phi=np.linspace(phi_os+pi,phi_ie-theta-2.0*pi*Nc,1000)
        (xi,yi)=coords_inv(phi, geo, theta, "fi")
        phi=np.linspace(phi_ie-theta-2.0*pi*Nc-pi,phi_os,1000)
        (xo,yo)=coords_inv(phi, geo, theta, "oo")
        V_poly=hs*polyarea(np.r_[xi,xo], np.r_[yi,yo])
        return Vd1,dVd1,V_poly
    
cpdef dict D1_forces(double theta, geoVals geo, bint poly = False):
    
    cython.declare(Nc = cython.long)
    
    hs=geo.h
    rb=geo.rb
    phi_ie=geo.phi_ie
    phi_e=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    phi_is=geo.phi_is
    phi_os=geo.phi_os
    ro=rb*(pi-phi_i0+phi_o0)
    Nc=getNc(theta, geo)
    
    #This is right before the discharge angle
    if abs(theta-theta_d(geo))<1e-8:
        return C1_forces(theta,Nc,geo,poly)
    
    phi2=phi_ie-theta-2.0*pi*Nc
    phi1=phi_os+pi
    VO=hs*rb**2/6.0*((phi2-phi_i0)**3-(phi1-phi_i0)**3)
    dVO=-hs*rb**2/2.0*((phi2-phi_i0)**2)
    
    phi2=phi_ie-theta-2.0*pi*Nc-pi
    phi1=phi_os
    VIa=hs*rb**2/6.0*((phi2-phi_o0)**3-(phi1-phi_o0)**3)
    dVIa=-hs*rb**2/2.0*((phi2-phi_o0)**2)
    
    VIb=hs*rb*ro/2.0*((phi_os-phi_o0)*sin(theta+phi_os-phi_ie)+cos(theta+phi_os-phi_ie))
    dVIb=hs*rb*ro/2.0*((phi_os-phi_o0)*cos(theta+phi_os-phi_ie)-sin(theta+phi_os-phi_ie))
    
    VIc=hs*rb*ro/2.0
    dVIc=0.0
    
    VId= hs*rb*ro/2.0*((phi_os-phi_i0+pi)*sin(theta+phi_os-phi_ie)+cos(theta+phi_os-phi_ie)+1)
    dVId=hs*rb*ro/2.0*((phi_os-phi_i0+pi)*cos(theta+phi_os-phi_ie)-sin(theta+phi_os-phi_ie))
    
    VI=VIa+VIb+VIc+VId
    dVI=dVIa+dVIb+dVIc+dVId
    
    Vd1=VO-VI
    dVd1=dVO-dVI

    cx_O=hs/VO*(fxA(rb,phi2,phi_i0)-fxA(rb,phi1,phi_i0))
    cy_O=hs/VO*(fyA(rb,phi2,phi_i0)-fyA(rb,phi1,phi_i0))
    cx_Ia=hs/VIa*(fxA(rb,phi2,phi_o0)-fxA(rb,phi1,phi_o0))
    cy_Ia=hs/VIa*(fyA(rb,phi2,phi_o0)-fyA(rb,phi1,phi_o0))
    cx_Ib=1.0/3.0*(-ro*sin(theta-phi_ie)+rb*(phi_os-phi_o0)*sin(phi_os)+rb*cos(phi_os))
    cy_Ib=1.0/3.0*(-ro*cos(theta-phi_ie)-rb*(phi_os-phi_o0)*cos(phi_os)+rb*sin(phi_os))
    cx_Ic=1.0/3.0*((rb*(-theta+phi_ie-phi_o0-2*pi*Nc-pi)-ro)*sin(theta-phi_ie)-rb*cos(theta-phi_ie))
    cy_Ic=1.0/3.0*((rb*(-theta+phi_ie-phi_o0-2*pi*Nc-pi)-ro)*cos(theta-phi_ie)+rb*sin(theta-phi_ie))
    cx_Id=(rb*(2*phi_os-phi_o0-phi_i0+pi)*sin(phi_os)-2*(ro*sin(theta-phi_ie)-rb*cos(phi_os)))/3.0
    cy_Id=(-2*(ro*cos(theta-phi_ie)-rb*sin(phi_os))-rb*(2*phi_os-phi_o0-phi_i0+pi)*cos(phi_os))/3.0
    cx_I=-(cx_Ia*VIa+cx_Ib*VIb+cx_Ic*VIc+cx_Id*VId)/VI+ro*cos(phi_ie-pi/2.0-theta)
    cy_I=-(cy_Ia*VIa+cy_Ib*VIb+cy_Ic*VIc+cy_Id*VId)/VI+ro*sin(phi_ie-pi/2.0-theta)
    
    cx=(cx_O*VO-cx_I*VI)/Vd1
    cy=(cy_O*VO-cy_I*VI)/Vd1
    
    fx_p=rb*hs*(sin(theta-phi_e)+(-theta-phi_o0+phi_e-2*pi*Nc-pi)*cos(theta-phi_e)-sin(phi_os)-(phi_o0-phi_os)*cos(phi_os))
    fy_p=-rb*hs*((-theta-phi_o0+phi_e-2*pi*Nc-pi)*sin(theta-phi_e)-cos(theta-phi_e)-(phi_os-phi_o0)*sin(phi_os)-cos(phi_os))
    M_O_p=(hs*rb**2*(theta-phi_os+2*phi_o0-phi_e+2*pi*Nc+pi)*(theta+phi_os-phi_e+2*pi*Nc+pi))/2.0
    fz_p = Vd1/hs
    
    exact_dict = dict(fx_p = fx_p,
                      fy_p = fy_p,
                      fz_p = fz_p,
                      M_O_p = M_O_p,
                      cx = cx,
                      cy = cy
                      )
    
    if not poly:
        return exact_dict
    else:
        ######################### Polygon calculations ##################
        phi=np.linspace(phi_os+pi,phi_ie-theta-2.0*pi*Nc,1000)
        (xi,yi)=coords_inv(phi, geo, theta, "fi")
        phi=np.linspace(phi_ie-theta-2.0*pi*Nc-pi,phi_os,1000)
        (xo,yo)=coords_inv(phi, geo, theta, "oo")
        V_poly=hs*polyarea(np.r_[xi,xo], np.r_[yi,yo])
        if V_poly>0:
            (cx_poly,cy_poly)=polycentroid(np.r_[xi,xo,xi[0]], np.r_[yi,yo,yi[0]])
        else:
            (cx_poly,cy_poly)=(xi[0], yi[0])
        ##################### Force Calculations #########################
        phi=np.linspace(phi_os,phi_ie-theta-2.0*pi*Nc-pi,1000)
        (xo,yo)=coords_inv(phi, geo, theta, "oo")
        nx=np.zeros_like(phi)
        ny=np.zeros_like(phi)
        (nx,ny)=coords_norm(phi,geo,theta,'oo')
        L=len(xo)
        dA=hs*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
        dfxp_poly=dA*(nx[1:L]+nx[0:L-1])/2.0
        dfyp_poly=dA*(ny[1:L]+ny[0:L-1])/2.0
        fxp_poly=np.sum(dfxp_poly)
        fyp_poly=np.sum(dfyp_poly)
        rOx=xo-geo.ro*cos(phi_e-pi/2-theta)
        rOx=(rOx[1:L]+rOx[0:L-1])/2
        rOy=yo-geo.ro*sin(phi_e-pi/2-theta)
        rOy=(rOy[1:L]+rOy[0:L-1])/2
        MO_poly=np.sum(rOx*dfyp_poly-rOy*dfxp_poly)
        poly_dict = dict(MO_poly = MO_poly,
                         fxp_poly = fxp_poly,
                         fyp_poly = fyp_poly,
                         cx_poly = cx_poly,
                         cy_poly = cy_poly
                         )
        exact_dict.update(poly_dict)
        return exact_dict

cpdef tuple D2(double theta, geoVals geo, bint poly=False):
    """
    Volume terms for D2 chamber
    
    Parameters
    ----------
    theta : float
        The crank angle in the range [:math:`0,2\pi`]
    geo : geoVals instance
        The geometry class
    poly : boolean, optional
        If true, also output the polygon calculations to the dict (SLOW!!)
    
    Returns
    -------
    values : tuple
        A tuple with volume,derivative of volume and volume from polygon(if requested)
    """
    #Is the same as the D1 chamber by symmetry
    return D1(theta,geo,poly)

cpdef dict D2_forces(double theta, geoVals geo, bint poly = False):

    cython.declare(Nc = cython.long)
    
    hs=geo.h
    rb=geo.rb
    phi_ie=geo.phi_ie
    phi_e=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    phi_is=geo.phi_is
    phi_os=geo.phi_os
    ro=rb*(pi-phi_i0+phi_o0)
    Nc=getNc(theta,geo=geo)
    
    #This is right before the discharge angle
    if abs(theta-theta_d(geo))<1e-8:
        return C2_forces(theta,Nc,geo,poly)

    phi2=phi_ie-theta-2.0*pi*Nc
    phi1=phi_os+pi
    VO=hs*rb**2/6.0*((phi2-phi_i0)**3-(phi1-phi_i0)**3)
    dVO=-hs*rb**2/2.0*((phi2-phi_i0)**2)
    
    phi2=phi_ie-theta-2.0*pi*Nc-pi
    phi1=phi_os
    VIa=hs*rb**2/6.0*((phi2-phi_o0)**3-(phi1-phi_o0)**3)
    dVIa=-hs*rb**2/2.0*((phi2-phi_o0)**2)
    
    VIb=hs*rb*ro/2.0*((phi_os-phi_o0)*sin(theta+phi_os-phi_ie)+cos(theta+phi_os-phi_ie))
    dVIb=hs*rb*ro/2.0*((phi_os-phi_o0)*cos(theta+phi_os-phi_ie)-sin(theta+phi_os-phi_ie))
    
    VIc=hs*rb*ro/2.0
    dVIc=0.0
    
    VId= hs*rb*ro/2.0*((phi_os-phi_i0+pi)*sin(theta+phi_os-phi_ie)+cos(theta+phi_os-phi_ie)+1)
    dVId=hs*rb*ro/2.0*((phi_os-phi_i0+pi)*cos(theta+phi_os-phi_ie)-sin(theta+phi_os-phi_ie))
    
    VI=VIa+VIb+VIc+VId
    dVI=dVIa+dVIb+dVIc+dVId
    
    Vd1=VO-VI
    dVd1=dVO-dVI

    cx_O=hs/VO*(fxA(rb,phi2,phi_i0)-fxA(rb,phi1,phi_i0))
    cy_O=hs/VO*(fyA(rb,phi2,phi_i0)-fyA(rb,phi1,phi_i0))
    cx_Ia=hs/VIa*(fxA(rb,phi2,phi_o0)-fxA(rb,phi1,phi_o0))
    cy_Ia=hs/VIa*(fyA(rb,phi2,phi_o0)-fyA(rb,phi1,phi_o0))
    cx_Ib=1.0/3.0*(-ro*sin(theta-phi_ie)+rb*(phi_os-phi_o0)*sin(phi_os)+rb*cos(phi_os))
    cy_Ib=1.0/3.0*(-ro*cos(theta-phi_ie)-rb*(phi_os-phi_o0)*cos(phi_os)+rb*sin(phi_os))
    cx_Ic=1.0/3.0*((rb*(-theta+phi_ie-phi_o0-2*pi*Nc-pi)-ro)*sin(theta-phi_ie)-rb*cos(theta-phi_ie))
    cy_Ic=1.0/3.0*((rb*(-theta+phi_ie-phi_o0-2*pi*Nc-pi)-ro)*cos(theta-phi_ie)+rb*sin(theta-phi_ie))
    cx_Id=(rb*(2*phi_os-phi_o0-phi_i0+pi)*sin(phi_os)-2*(ro*sin(theta-phi_ie)-rb*cos(phi_os)))/3.0
    cy_Id=(-2*(ro*cos(theta-phi_ie)-rb*sin(phi_os))-rb*(2*phi_os-phi_o0-phi_i0+pi)*cos(phi_os))/3.0
    cx_I=-(cx_Ia*VIa+cx_Ib*VIb+cx_Ic*VIc+cx_Id*VId)/VI+ro*cos(phi_ie-pi/2.0-theta)
    cy_I=-(cy_Ia*VIa+cy_Ib*VIb+cy_Ic*VIc+cy_Id*VId)/VI+ro*sin(phi_ie-pi/2.0-theta)
    
    cxd1=(cx_O*VO-cx_I*VI)/Vd1
    cyd1=(cy_O*VO-cy_I*VI)/Vd1
    
    fx_p=-hs*rb*(-sin(theta-phi_e)+(theta+phi_i0-phi_e+2*pi*Nc)*cos(theta-phi_e)+sin(phi_os)-(phi_os-phi_i0+pi)*cos(phi_os))
    fy_p=hs*rb*((theta+phi_i0-phi_e+2*pi*Nc)*sin(theta-phi_e)+cos(theta-phi_e)-(-phi_os+phi_i0-pi)*sin(phi_os)+cos(phi_os))
    M_O_p=-(hs*rb**2*(theta-phi_os+2*phi_i0-phi_e+2*pi*Nc-pi)*(theta+phi_os-phi_e+2*pi*Nc+pi))/2
    
    (cx,cy)=(-cxd1+ro*cos(phi_ie-pi/2-theta),-cyd1+ro*sin(phi_ie-pi/2-theta))
    fz_p = Vd1/hs
    
    exact_dict = dict(fx_p = fx_p,
                      fy_p = fy_p,
                      fz_p = fz_p,
                      M_O_p = M_O_p,
                      cx = cx,
                      cy = cy
                      )
    
    if not poly:
        return exact_dict
    else:
        raise NotImplementedError('D2_forces polygon not implemented')
#        phi=np.linspace(phi_os+pi,phi_ie-theta-2.0*pi*Nc,1000)
#        (xo,yo)=coords_inv(phi, geo, theta, "oi")
#        nx=np.zeros_like(phi)
#        ny=np.zeros_like(phi)
#        (nx,ny)=coords_norm(phi,geo,theta,"oi")
#        L=len(xo)
#        dA=h*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
#        fxp_poly=np.sum(dA*(nx[1:L]+nx[0:L-1])/2.0)
#        fyp_poly=np.sum(dA*(ny[1:L]+ny[0:L-1])/2.0)
#        return exact_dict
    
cpdef tuple DD(double theta, geoVals geo, bint poly=False):
    
    hs=geo.h
    xa1=geo.xa_arc1
    ya1=geo.ya_arc1
    ra1=geo.ra_arc1
    ta1_2=geo.t2_arc1
    ta1_1=geo.t1_arc1
    xa2=geo.xa_arc2
    ya2=geo.ya_arc2
    ra2=geo.ra_arc2
    ta2_2=geo.t2_arc2
    ta2_1=geo.t1_arc2    
    ro=geo.ro
    m_line=geo.m_line
    b_line=geo.b_line
    t1_line=geo.t1_line
    t2_line=geo.t2_line
    phi_os=geo.phi_os
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    phi_is=geo.phi_is
    phi_e=geo.phi_ie
    rb=geo.rb
    om=geo.phi_ie-pi/2-theta
        
    (xoos,yoos)=coords_inv(geo.phi_os, geo, theta, 'oo')
    
    #################### Oa portion ####################        
    V_Oa=hs*((-(ra1*(cos(ta1_2)*(ya1-yoos)-sin(ta1_2)*(xa1-xoos)-ra1*ta1_2))/2)-(-(ra1*(cos(ta1_1)*(ya1-yoos)-sin(ta1_1)*(xa1-xoos)-ra1*ta1_1))/2))
    dV_Oa=-hs*ra1*ro/2.0*((sin(om)*sin(ta1_2)+cos(om)*cos(ta1_2))-(sin(om)*sin(ta1_1)+cos(om)*cos(ta1_1)))
    
    #################### Ob portion ####################
    
    x1l=t1_line #old nomenclature
    y1l=m_line*t1_line+b_line #old nomenclature
    V_Ob=hs/2.0*((ro*xoos-ro*x1l)*sin(om)-(ro*cos(om)-2.0*x1l)*yoos+y1l*(ro*cos(om)-2.0*xoos))
    dV_Ob=ro*hs/2.0*(ro-yoos*sin(om)-xoos*cos(om)-y1l*sin(om)-x1l*cos(om))
    
    ##################### Oc portion ###################
    V_Oc=rb*hs/6*(
          3*ro*(phi_os-phi_i0+pi)*sin(theta+phi_os-phi_e)
          +3*ro*cos(theta+phi_os-phi_e)
          +3*(phi_is-phi_i0)*ro*sin(theta+phi_is-phi_e)
          +3*ro*cos(theta+phi_is-phi_e)
          +3*rb*((phi_is-phi_i0)*(phi_os-phi_o0)+1)*sin(phi_os-phi_is)
          -3*rb*(phi_os-phi_o0-phi_is+phi_i0)*cos(phi_os-phi_is)
          +rb*((phi_os+pi-phi_i0)**3-(phi_is-phi_i0)**3)+3*ro)
    dV_Oc=rb*hs*ro/2*(
           (phi_os-phi_i0+pi)*cos(theta+phi_os-phi_e)
          -sin(theta+phi_os-phi_e)
          +(phi_is-phi_i0)*cos(theta+phi_is-phi_e)
          -sin(theta+phi_is-phi_e)
          )

    #################### Ia portion ####################

    V_Ia=hs*ra2/2.0*(xa2*(sin(ta2_2)-sin(ta2_1))
                 -ya2*(cos(ta2_2)-cos(ta2_1))
                 -rb*(sin(ta2_2-phi_os)-sin(ta2_1-phi_os))
        -rb*(phi_os-phi_o0)*(cos(ta2_2-phi_os)-cos(ta2_1-phi_os))
                 +ra2*(ta2_2-ta2_1)  )
    dV_Ia=0.0
    
    #################### Ib portion #####################
    x1l=t1_line #old nomenclature
    x2l=t2_line #old nomenclature
    y1l=m_line*t1_line+b_line #old nomenclature
    ml=m_line
    V_Ib=-hs*(x2l-x1l)/2.0*(rb*ml*(cos(phi_os)+(phi_os-phi_o0)*sin(phi_os))+b_line-rb*(sin(phi_os)-(phi_os-phi_o0)*cos(phi_os)))
    dV_Ib=0
    
    cx=ro*cos(om)/2.0 #By symmetry
    cy=ro*sin(om)/2.0 #By symmetry
    V=2.0*(V_Oa+V_Ob+V_Oc-V_Ia-V_Ib)
    dV=2.0*(dV_Oa+dV_Ob+dV_Oc-dV_Ia-dV_Ib)
    
    if not poly:
        return V,dV
    else:
        ##########################################################
        ##                    POLYGON                           ##
        ##########################################################
        t=np.linspace(geo.t1_arc1,geo.t2_arc1,300)
        (x_farc1,y_farc1)=(
            geo.xa_arc1+geo.ra_arc1*cos(t),
            geo.ya_arc1+geo.ra_arc1*sin(t))
        (x_oarc1,y_oarc1)=(
           -geo.xa_arc1-geo.ra_arc1*cos(t)+geo.ro*cos(om),
           -geo.ya_arc1-geo.ra_arc1*sin(t)+geo.ro*sin(om))
        
        t=np.linspace(geo.t1_arc2,geo.t2_arc2,300)
        (x_farc2,y_farc2)=(
            geo.xa_arc2+geo.ra_arc2*np.cos(t),
            geo.ya_arc2+geo.ra_arc2*np.sin(t))
        (x_oarc2,y_oarc2)=(
           -geo.xa_arc2-geo.ra_arc2*np.cos(t)+geo.ro*cos(om),
           -geo.ya_arc2-geo.ra_arc2*np.sin(t)+geo.ro*sin(om)) 
        
        phi=np.linspace(phi_is,phi_os+pi,300)
        (x_finv,y_finv)=coords_inv(phi,geo,theta,'fi')
        (x_oinv,y_oinv)=coords_inv(phi,geo,theta,'oi')
        
        x=np.r_[x_farc2[::-1],x_farc1,x_finv,x_oarc2[::-1],x_oarc1,x_oinv,x_farc2[-1]]
        y=np.r_[y_farc2[::-1],y_farc1,y_finv,y_oarc2[::-1],y_oarc1,y_oinv,y_farc2[-1]]
        V_poly=geo.h*polyarea(x, y)
        return V,dV,V_poly
    
cpdef dict DD_forces(double theta, geoVals geo, bint poly=False):
        
    hs=geo.h
    xa1=geo.xa_arc1
    ya1=geo.ya_arc1
    ra1=geo.ra_arc1
    ta1_2=geo.t2_arc1
    ta1_1=geo.t1_arc1
    xa2=geo.xa_arc2
    ya2=geo.ya_arc2
    ra2=geo.ra_arc2
    ta2_2=geo.t2_arc2
    ta2_1=geo.t1_arc2    
    ro=geo.ro
    m_line=geo.m_line
    b_line=geo.b_line
    t1_line=geo.t1_line
    t2_line=geo.t2_line
    phi_os=geo.phi_os
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    phi_is=geo.phi_is
    phi_e=geo.phi_ie
    rb=geo.rb
    om=geo.phi_ie-pi/2-theta
    
    ################ Force Components #########
    #Arc 1
    fx_p =-hs*geo.ra_arc1*(sin(geo.t2_arc1)-sin(geo.t1_arc1))
    fy_p =+hs*geo.ra_arc1*(cos(geo.t2_arc1)-cos(geo.t1_arc1))
    M_O_p  =-hs*geo.ra_arc1*((sin(geo.t2_arc1)-sin(geo.t1_arc1))*geo.ya_arc1+(cos(geo.t2_arc1)-cos(geo.t1_arc1))*geo.xa_arc1)
    #Arc 2
    fx_p+=+hs*geo.ra_arc2*(sin(geo.t2_arc2)-sin(geo.t1_arc2))
    fy_p+=-hs*geo.ra_arc2*(cos(geo.t2_arc2)-cos(geo.t1_arc2))
    M_O_p +=+hs*geo.ra_arc2*((sin(geo.t2_arc2)-sin(geo.t1_arc2))*geo.ya_arc2+(cos(geo.t2_arc2)-cos(geo.t1_arc2))*geo.xa_arc2)
    
    #Line
    x1t=-geo.xa_arc1-geo.ra_arc1*cos(geo.t1_arc1)+ro*cos(om)
    y1t=-geo.ya_arc1-geo.ra_arc1*sin(geo.t1_arc1)+ro*sin(om)
    x2t=-geo.xa_arc2-geo.ra_arc2*cos(geo.t1_arc2)+ro*cos(om)
    y2t=-geo.ya_arc2-geo.ra_arc2*sin(geo.t1_arc2)+ro*sin(om)
    L=np.sqrt((x2t-x1t)**2+(y2t-y1t)**2)
    if L>1e-12:
        Lx=(x2t-x1t)/L
        Ly=(y2t-y1t)/L
        nx=-1/np.sqrt(1+Lx**2/Ly**2)
        ny=Lx/Ly/np.sqrt(1+Lx**2/Ly**2)
        # Make sure you get the cross product with the normal 
        # pointing towards the scroll, otherwise flip...
        if Lx*ny-Ly*nx<0:
            nx*=-1
            ny*=-1
        fx_p+=hs*nx*L
        fy_p+=hs*ny*L
        rx=(x1t+x2t)/2-ro*cos(om)
        ry=(y2t+y2t)/2-ro*sin(om)
        M_O_p+=rx*hs*ny*L-ry*hs*nx*L
    
    #Involute portion
    fx_p+=-hs*(-sin(phi_os)+(phi_os-phi_i0+pi)*cos(phi_os)-sin(phi_is)-(phi_i0-phi_is)*cos(phi_is))*rb
    fy_p+=hs*((-phi_os+phi_i0-pi)*sin(phi_os)-cos(phi_os)-(phi_is-phi_i0)*sin(phi_is)-cos(phi_is))*rb
    M_O_p +=-(hs*(phi_os-phi_is+pi)*(phi_os+phi_is-2*phi_i0+pi)*rb*rb)/2
    
    cx=ro*cos(om)/2.0
    cy=ro*sin(om)/2.0
    
    #Axial force
    fz_p = DD(theta,geo)[0]/hs
    
    exact_dict = dict(fx_p = fx_p,
                      fy_p = fy_p,
                      fz_p = fz_p,
                      M_O_p = M_O_p,
                      cx = cx,
                      cy = cy
                      )
    
    if not poly:
        return exact_dict
    else:
        ##########################################################
        ##                    POLYGON                           ##
        ##########################################################
        t=np.linspace(geo.t1_arc1,geo.t2_arc1,300)
        (x_farc1,y_farc1)=(
            geo.xa_arc1+geo.ra_arc1*np.cos(t),
            geo.ya_arc1+geo.ra_arc1*np.sin(t))
        (x_oarc1,y_oarc1)=(
           -geo.xa_arc1-geo.ra_arc1*np.cos(t)+geo.ro*cos(om),
           -geo.ya_arc1-geo.ra_arc1*np.sin(t)+geo.ro*sin(om))
        (nx_oarc1,ny_oarc1)=(-np.cos(t),-np.sin(t))
        
        t=np.linspace(geo.t1_arc2,geo.t2_arc2,300)
        (x_farc2,y_farc2)=(
            geo.xa_arc2+geo.ra_arc2*np.cos(t),
            geo.ya_arc2+geo.ra_arc2*np.sin(t))
        (x_oarc2,y_oarc2)=(
           -geo.xa_arc2-geo.ra_arc2*np.cos(t)+geo.ro*cos(om),
           -geo.ya_arc2-geo.ra_arc2*np.sin(t)+geo.ro*sin(om)) 
        (nx_oarc2,ny_oarc2)=(+np.cos(t),+np.sin(t))
        
        phi=np.linspace(phi_is,phi_os+pi,300)
        (x_finv,y_finv)=coords_inv(phi,geo,theta,'fi')
        (x_oinv,y_oinv)=coords_inv(phi,geo,theta,'oi')
        (nx_oinv,ny_oinv)=coords_norm(phi,geo,theta,'oi')
        
        x=np.r_[x_farc2[::-1],x_farc1,x_finv,x_oarc2[::-1],x_oarc1,x_oinv,x_farc2[-1]]
        y=np.r_[y_farc2[::-1],y_farc1,y_finv,y_oarc2[::-1],y_oarc1,y_oinv,y_farc2[-1]]
        (cx_poly,cy_poly)=polycentroid(x,y)
        V_poly=geo.h*polyarea(x, y)
        
        fxp_poly=0
        fyp_poly=0
        MO_poly=0
        #Arc1
        L=len(nx_oarc1)
        dA=hs*np.sqrt(np.power(x_oarc1[1:L]-x_oarc1[0:L-1],2)+np.power(y_oarc1[1:L]-y_oarc1[0:L-1],2))
        dfxp_poly=dA*(nx_oarc1[1:L]+nx_oarc1[0:L-1])/2.0
        dfyp_poly=dA*(ny_oarc1[1:L]+ny_oarc1[0:L-1])/2.0
        fxp_poly=np.sum(dfxp_poly)
        fyp_poly=np.sum(dfyp_poly)
        rOx=x_oarc1-geo.ro*cos(phi_e-pi/2-theta)
        rOx=(rOx[1:L]+rOx[0:L-1])/2
        rOy=y_oarc1-geo.ro*sin(phi_e-pi/2-theta)
        rOy=(rOy[1:L]+rOy[0:L-1])/2
        MO_poly=np.sum(rOx*dfyp_poly-rOy*dfxp_poly)
        print 'Arc1',np.sum(dfxp_poly),np.sum(dfyp_poly)
        #Arc2
        L=len(nx_oarc2)
        dA=hs*np.sqrt(np.power(x_oarc2[1:L]-x_oarc2[0:L-1],2)+np.power(y_oarc2[1:L]-y_oarc2[0:L-1],2))
        dfxp_poly=dA*(nx_oarc2[1:L]+nx_oarc2[0:L-1])/2.0
        dfyp_poly=dA*(ny_oarc2[1:L]+ny_oarc2[0:L-1])/2.0
        fxp_poly+=np.sum(dfxp_poly)
        fyp_poly+=np.sum(dfyp_poly)
        rOx=x_oarc2-geo.ro*cos(phi_e-pi/2-theta)
        rOx=(rOx[1:L]+rOx[0:L-1])/2
        rOy=y_oarc2-geo.ro*sin(phi_e-pi/2-theta)
        rOy=(rOy[1:L]+rOy[0:L-1])/2
        MO_poly+=np.sum(rOx*dfyp_poly-rOy*dfxp_poly)
        print 'Arc2',np.sum(dfxp_poly),np.sum(dfyp_poly)
        #Involute
        L=len(y_oinv)
        dA=hs*np.sqrt(np.power(x_oinv[1:L]-x_oinv[0:L-1],2)+np.power(y_oinv[1:L]-y_oinv[0:L-1],2))
        dfxp_poly=dA*(nx_oinv[1:L]+nx_oinv[0:L-1])/2.0
        dfyp_poly=dA*(ny_oinv[1:L]+ny_oinv[0:L-1])/2.0
        fxp_poly+=np.sum(dfxp_poly)
        fyp_poly+=np.sum(dfyp_poly)
        rOx=x_oinv-geo.ro*cos(phi_e-pi/2-theta)
        rOx=(rOx[1:L]+rOx[0:L-1])/2
        rOy=y_oinv-geo.ro*sin(phi_e-pi/2-theta)
        rOy=(rOy[1:L]+rOy[0:L-1])/2
        MO_poly+=np.sum(rOx*dfyp_poly-rOy*dfxp_poly)
        print 'Involute',np.sum(dfxp_poly),np.sum(dfyp_poly)
        #Line
        x1t=-geo.xa_arc1-geo.ra_arc1*cos(geo.t1_arc1)+ro*cos(om)
        y1t=-geo.ya_arc1-geo.ra_arc1*sin(geo.t1_arc1)+ro*sin(om)
        x2t=-geo.xa_arc2-geo.ra_arc2*cos(geo.t1_arc2)+ro*cos(om)
        y2t=-geo.ya_arc2-geo.ra_arc2*sin(geo.t1_arc2)+ro*sin(om)
        L=np.sqrt((x2t-x1t)**2+(y2t-y1t)**2)
        if L>1e-12:
            Lx=(x2t-x1t)/L
            Ly=(y2t-y1t)/L
            nx=-1/np.sqrt(1+Lx**2/Ly**2)
            ny=Lx/Ly/np.sqrt(1+Lx**2/Ly**2)
            # Make sure you get the cross product with the normal 
            # pointing towards the scroll, otherwise flip...
            if Lx*ny-Ly*nx<0:
                nx*=-1
                ny*=-1
            fxp_poly+=hs*nx*L
            fyp_poly+=hs*ny*L
            rx=(x1t+x2t)/2-ro*cos(om)
            ry=(y2t+y2t)/2-ro*sin(om)
            MO_poly+=rx*hs*ny*L-ry*hs*nx*L
            print 'Line',hs*nx*L,hs*ny*L
            
        poly_dict = dict(MO_poly = MO_poly,
                         fxp_poly = fxp_poly,
                         fyp_poly = fyp_poly,
                         cx_poly = cx_poly,
                         cy_poly = cy_poly
                         )
        exact_dict.update(poly_dict)
        return exact_dict
   
cpdef tuple DDD(double theta, geoVals geo, bint poly=False): 
    
    if not poly:
        V_d1,dV_d1=D1(theta,geo)
        V_d2,dV_d2=D2(theta,geo)
        V_dd,dV_dd=DD(theta,geo)
        V_ddd =   V_d1+ V_d2+ V_dd
        dV_ddd = dV_d1+dV_d2+dV_dd
        return V_ddd,dV_ddd
    else:
        raise AttributeError('Polygons not coded for DDD chamber')

cpdef dict DDD_forces(double theta, geoVals geo, bint poly=False):
    
    if not poly:
        exact_dict = {}
        _D1_forces = D1_forces(theta,geo)
        _D2_forces = D2_forces(theta,geo)
        _DD_forces = DD_forces(theta,geo)
        for key in _D1_forces:
            exact_dict[key] = _D1_forces[key]+_D2_forces[key]+_DD_forces[key]
        return exact_dict
    else:
        raise AttributeError('Polygons not coded for DDD chamber')
        
 
cpdef double phi_s_sa(double theta, geoVals geo):
    
    h=geo.h
    rb=geo.rb 
    phi_ie=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    ro=rb*(pi-phi_i0+phi_o0)
    
    b=(-phi_o0+phi_ie-pi)
    D=ro/rb*((phi_i0-phi_ie)*sin(theta)-cos(theta)+1)/(phi_ie-phi_i0)
    B=1.0/2.0*(sqrt(b**2-4.0*D)-b)
    return phi_ie-pi+B-phi_o0
    
cpdef double phi_d_dd(double theta, geoVals geo):
    
    phi_os=geo.phi_os;
    phi_o0=geo.phi_o0;
    phi_ie=geo.phi_ie;
    phi_i0=geo.phi_i0;
    alpha=pi-phi_i0+phi_o0;
    
    # Use secant method to calculate the involute angle at break
    eps=1e-8;
    change=999;
    iter=1;
    while ((iter<=3 or abs(f)>eps) and iter<100):
        if (iter==1): x1=geo.phi_is; phi=x1;
        if (iter==2): x2=geo.phi_is+0.1; phi=x2;
        if (iter>2): phi=x2;

        f=1+cos(phi-phi_os)-(phi_os-phi_o0)*sin(phi-phi_os)+alpha*sin(phi-phi_ie+theta);

        if (iter==1): y1=f;
        if (iter==2): y2=f;
        if (iter>2):
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
            
        iter+=1
        
        # If the value is still less than the starting angle
        # after 20 iterations
        if (iter>20 and x3<geo.phi_is):
            return geo.phi_is;

    if (x3>geo.phi_is): return x3;
    else: return geo.phi_is;

cpdef double Area_d_dd(double theta, geoVals geo):
    x_fis,y_fis=coords_inv(phi_d_dd(theta,geo),geo,theta,"fi")
    x_oos,y_oos=coords_inv(geo.phi_os,geo,theta,"oo")
    return geo.h*((x_fis-x_oos)**2+(y_fis-y_oos)**2)**0.5

cpdef double Area_s_s1_offset(double theta, geoVals geo):
    """
    The area between the suction area and the s1 chamber when an offset scroll
    wrap is employed [:math:`m^2`]
    
    Parameters
    ----------
    theta : float
    geo : geoVals instance
    
    Returns
    -------
    A : float
        Area [:math:`m^2`]
    """
    if theta < pi:
        w = geo.ro*(1+cos(theta)) + geo.delta_suction_offset
    else:
        w = geo.delta_suction_offset
    return geo.h * w
            
cpdef double Area_s_sa(double theta, geoVals geo):
    """
    The area between the suction area and the s1 chamber when an offset scroll
    wrap is not employed and the area between the suction area and the 
    s2 chamber [:math:`m^2`]
    
    Parameters
    ----------
    theta : float
    geo : geoVals instance
    
    Returns
    -------
    A : float
        Area [:math:`m^2`]
    """
    return geo.ro*geo.h*(1-cos(theta))
    
if __name__=='__main__':
    """
    xo : ro*cos(phi_ie-%pi/2-%theta)$
    yo : ro*sin(phi_ie-%pi/2-%theta)$
    xssa : rb*(cos(phi_ie-%pi+B+delta_phi_ie)+(phi_ie-%pi+B+delta_phi_ie-phi_o0)*sin(phi_ie-%pi+B+delta_phi_ie))$
    yssa : rb*(sin(phi_ie-%pi+B+delta_phi_ie)-(phi_ie-%pi+B+delta_phi_ie-phi_o0)*cos(phi_ie-%pi+B+delta_phi_ie))$
    A : -(xo*yssa-yo*xssa)/2$
    trigsimp(%)$
    facsum(%,[sin(B+phi_ie),cos(B+phi_ie),sin(phi_ie-%theta),
        cos(%theta-phi_ie),[phi_ie,phi_0,theta]])$
    trigsimp(trigreduce(expand(%)));
    """
    print 'This is the base file with scroll geometry.  Running this file doesn\'t do anything'