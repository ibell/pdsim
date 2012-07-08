
from __future__ import division

import numpy as np
import cython
import time #(temporary)
 
#This is a list of all the members in geoVals
geoValsvarlist=['h','phi_i0','phi_is','phi_ie','phi_e','phi_o0','ro','rb','phi_os','phi_oe','t',
'xa_arc1','ya_arc1','ra_arc1','t1_arc1','t2_arc1',
'xa_arc2','ya_arc2','ra_arc2','t1_arc2','t2_arc2',
'b_line', 't1_line', 't2_line', 'm_line',
'x0_wall','y0_wall','r_wall',
'delta_radial', 'delta_flank']
 
def rebuild_geoVals(d):
    geo = geoVals()
    for atr in geoValsvarlist:
        setattr(geo,atr,d[atr])
    return geo 
    
class geoVals:
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
        
def fxA(rb,phi,phi0):
    return rb**3/3.0*(4.0*((phi-phi0)**2-2.0)*sin(phi)+(phi0-phi)*((phi-phi0)**2-8.0)*cos(phi))

def fyA(rb,phi,phi0):
    return rb**3/3.0*((phi0-phi)*((phi-phi0)**2-8.0)*sin(phi)-4.0*((phi-phi0)**2-2.0)*cos(phi))

def theta_d(geo):  
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

def nC_Max(geo):
    return int(floor((geo.phi_ie-geo.phi_os-pi)/(2.0*pi)))

def getNc(theta,geo):
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
    
def _coords_inv_np(phi,geo,theta,flag=""):
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
    
def _coords_inv_d(phi,geo,theta,flag=""):
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
  
def coords_inv(phi,geo,theta,flag="fi"):
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
            raise AttributeError('error with starting angles phi_os %.16f phi_is-pi %.16f' %(geo.phi_os,geo.phi_is-pi))
            
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
  
def SA(theta, geo, poly=False, forces=False):
    h=geo.h
    rb=geo.rb 
    phi_ie=geo.phi_ie 
    phi_o0=geo.phi_o0
    phi_oe=geo.phi_oe
    phi_i0=geo.phi_i0
    oe=geo.phi_oe
    ro=rb*(pi-phi_i0+phi_o0)
    
    b=(-phi_o0+phi_ie-pi);
    D=ro/rb*((phi_i0-phi_ie)*sin(theta)-cos(theta)+1)/(phi_ie-phi_i0);
    B=1.0/2.0*(sqrt(b*b-4.0*D)-b);
    B_prime=-ro/rb*(sin(theta)+(phi_i0-phi_ie)*cos(theta))/((phi_ie-phi_i0)*sqrt(b*b-4*D));
    V_Isa=h*rb**2/6.0*(pow(phi_oe-phi_o0,3)-pow(phi_ie-pi+B-phi_o0,3));
    
    V=h*pi*geo.r_wall**2-2*V_Isa;
    dV=h*rb**2*pow(phi_ie-pi+B-phi_o0,2)*B_prime;
    
    if forces==False and poly==False:
        return V,dV
    else:
        raise ValueError('no forces or polygons for SA yet')
    
def S1(theta, geo, poly=False, forces=False, theta_0_volume=1e-9):
        
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
    
    VIc=h*rb*ro/2
    dVIc=0
        
    Vs=VO-(VIa+VIb-VIc)
    dVs=dVO-(dVIa+dVIb-dVIc)
    
    #Add on the ficticious volume to correct for there actually being no volume
    #  at theta=0
    Vs+=theta_0_volume
    
    if forces==False and poly==False:
        return Vs,dVs
        
    if forces==True:
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
        
        fx_p=-rb*h*(sin(B+phi_e)-(B-phi_o0+phi_e-pi)*cos(B+phi_e)+sin(theta-phi_e)-(theta+phi_o0-phi_e+pi)*cos(theta-phi_e))
        fy_p=rb*h*((B-phi_o0+phi_e-pi)*sin(B+phi_e)+cos(B+phi_e)-(theta+phi_o0-phi_e+pi)*sin(theta-phi_e)-cos(theta-phi_e))
        M_O=(h*rb**2*(B-theta-2*phi_o0+2*phi_e-2*pi)*(B+theta))/2
    else:
        (cx,cy,fx_p,fy_p,M_O)=(None,None,None,None,None)
    
    if poly==True:
        ############### Polygon calculations ##################
        phi=np.linspace(phi_ie-theta,phi_ie,2000)
        (xi,yi)=coords_inv(phi, geo, theta, 'fi')
        phi=np.linspace(phi_ie-pi+B,phi_ie-pi-theta,2000)
        (xo,yo)=coords_inv(phi, geo, theta, 'oo')
        V_poly=h*polyarea(np.r_[xi,xo,xi[0]], np.r_[yi,yo,yi[0]])
        (cx_poly,cy_poly)=polycentroid(np.r_[xi,xo,xi[0]], np.r_[yi,yo,yi[0]])
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
    else:
        (V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly)=(None,None,None,None,None,None)
    
    if forces==False and poly==True:
        return Vs,dVs,V_poly
    else:
        return Vs,dVs,cx,cy,fx_p,fy_p,M_O,V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly

def S2(theta, geo,poly=False, forces=False, theta_0_volume=1e-9):
    
    if forces==False and poly==False:
        return S1(theta,geo,theta_0_volume=theta_0_volume)
    else:
        h=geo.h
        rb=geo.rb
        phi_ie=geo.phi_ie
        phi_e=geo.phi_ie
        phi_o0=geo.phi_o0
        phi_i0=geo.phi_i0
        ro=rb*(pi-phi_i0+phi_o0)
        (Vs1,dVs1,cx_s1,cy_s1,fx_ps1,fy_ps1,M_O_s1,cs1,V_polys1,cx_polys1,cy_polys1,fxp_polys1,fyp_polys1,M_Os1_poly)=S1(theta,geo)
        (cx,cy)=(-cx_s1+ro*cos(phi_ie-pi/2-theta),-cy_s1+ro*sin(phi_ie-pi/2-theta))
    
        fx_p=-rb*h*(sin(theta-phi_e)-(theta+phi_i0-phi_e)*cos(theta-phi_e)+cos(phi_e)*(phi_i0-phi_e)+sin(phi_e))
        fy_p=-rb*h*((theta+phi_i0-phi_e)*sin(theta-phi_e)+cos(theta-phi_e)+sin(phi_e)*(phi_i0-phi_e)-cos(phi_e))
        M_O=(h*rb**2*theta*(theta+2*phi_i0-2*phi_e))/2
    
    if poly==True:
        ############### Numerical Force Calculations ###########
        phi=np.linspace(phi_ie-theta,phi_ie,2000)
        (xo,yo)=coords_inv(phi, geo, theta, 'oi')
        nx=np.zeros_like(phi)
        ny=np.zeros_like(phi)
        (nx,ny)=coords_norm(phi,geo,theta,'oi')
        L=len(xo)
        dA=h*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
        fxp_poly=np.sum(dA*(nx[1:L]+nx[0:L-1])/2.0)
        fyp_poly=np.sum(dA*(ny[1:L]+ny[0:L-1])/2.0)
    else:
        (V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly)=(None,None,None,None,None,None)
    
    return (Vs1,dVs1,cx,cy,fx_p,fy_p,fxp_poly,fyp_poly)

def C1(theta, alpha, geo, poly=False, forces=False):
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
    if forces==True:
        psi=rb/3.0*(3.0*theta**2+6.0*phi_o0*theta+3.0*phi_o0**2+pi**2-15.0+(theta+phi_o0)*(12.0*pi*alpha-6.0*phi_ie)+3.0*phi_ie**2+12.0*pi*alpha*(pi*alpha-phi_ie))/(2.0*theta+phi_o0-2.0*phi_ie+phi_i0+4.0*pi*alpha-pi)
        cx=-2.0*rb*cos(theta-phi_ie)-psi*sin(theta-phi_ie)
        cy=+2.0*rb*sin(theta-phi_ie)-psi*cos(theta-phi_ie) 
        fx_p= 2.0*pi*rb*h*cos(theta-phi_e)
        fy_p=-2.0*pi*rb*h*sin(theta-phi_e)
        M_O=-2*pi*h*rb*rb*(theta+phi_o0-phi_e+2*pi*alpha)
    
    if poly==True:
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
    else:
        (V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly)=(None,None,None,None,None,None)
    
    if forces==False and poly==False:
        return V,dV
    elif forces==False and poly==True:
        return V,dV,V_poly
    else:
        return (V,dV,cx,cy,fx_p,fy_p,M_O,V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly)

def C2(theta, alpha, geo, poly=False, forces=False):
    
    if forces==False and poly==False:
        return C1(theta,alpha,geo)
    else:
        h=geo.h
        rb=geo.rb
        phi_ie=geo.phi_ie
        phi_e=geo.phi_ie
        phi_o0=geo.phi_o0
        phi_i0=geo.phi_i0
        ro=rb*(pi-phi_i0+phi_o0)
        
        ro=geo.ro
        phi_ie=geo.phi_ie
        (Vc1,dVc1,cxc1,cyc1,fx_pc1,fy_pc1,M_Oc1,V_polyc1,cx_polyc1,cy_polyc1,fxp_polyc1,fyp_polyc1,M_Oc1_poly)=C1(theta,alpha,geo)
        (cx,cy)=(-cxc1+ro*cos(phi_ie-pi/2-theta),-cyc1+ro*sin(phi_ie-pi/2-theta))
        fx_p= 2.0*pi*rb*h*cos(theta-phi_e)
        fy_p=-2.0*pi*rb*h*sin(theta-phi_e)
        M_O=2*pi*h*rb*rb*(theta+phi_i0-phi_e+2*pi*alpha-pi)
    
    if poly==True:
        ##################### Force Calculations #########################
        phi=np.linspace( geo.phi_ie-theta-2*pi*(alpha),geo.phi_ie-theta-2*pi*(alpha-1),1000)
        (xo,yo)=coords_inv(phi, geo, theta, 'oi')
        nx=np.zeros_like(phi)
        ny=np.zeros_like(phi)
        (nx,ny)=coords_norm(phi,geo,theta,'oi')
        L=len(xo)
        dA=h*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
        fxp_poly=np.sum(dA*(nx[1:L]+nx[0:L-1])/2.0)
        fyp_poly=np.sum(dA*(ny[1:L]+ny[0:L-1])/2.0)
    else:
        (V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly)=(None,None,None,None,None,None)
        
    return (Vc1,dVc1,cx,cy,fx_p,fy_p,fxp_poly,fyp_poly)

def D1(theta, geo, poly=False, forces=False):
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
    dVIc=0
    
    VId= hs*rb*ro/2.0*((phi_os-phi_i0+pi)*sin(theta+phi_os-phi_ie)+cos(theta+phi_os-phi_ie)+1)
    dVId=hs*rb*ro/2.0*((phi_os-phi_i0+pi)*cos(theta+phi_os-phi_ie)-sin(theta+phi_os-phi_ie))
    
    VI=VIa+VIb+VIc+VId
    dVI=dVIa+dVIb+dVIc+dVId
    
    Vd1=VO-VI
    dVd1=dVO-dVI
    
    if forces==False and poly==False:
        return Vd1,dVd1
    
    if forces==True:
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
        M_O=(hs*rb**2*(theta-phi_os+2*phi_o0-phi_e+2*pi*Nc+pi)*(theta+phi_os-phi_e+2*pi*Nc+pi))/2.0
    
    if poly==True:
        ######################### Polygon calculations ##################
        phi=np.linspace(phi_os+pi,phi_ie-theta-2.0*pi*Nc,1000)
        (xi,yi)=coords_inv(phi, geo, theta, "fi")
        phi=np.linspace(phi_ie-theta-2.0*pi*Nc-pi,phi_os,1000)
        (xo,yo)=coords_inv(phi, geo, theta, "oo")
        V_poly=hs*polyarea(np.r_[xi,xo], np.r_[yi,yo])
        (cx_poly,cy_poly)=polycentroid(np.r_[xi,xo,xi[0]], np.r_[yi,yo,yi[0]])
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
    else:
        (V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly)=(None,None,None,None,None,None)
         
    return Vd1,dVd1,cx,cy,fx_p,fy_p,M_O,(cx_Ia,cy_Ia,cx_Ib,cy_Ib,cx_Ic,cy_Ic,cx_Id,cy_Id,cx_I,cy_I,cx_O,cy_O,VO,VIa,VIb,VIc,VId,dVO,dVIa,dVIb,dVIc,dVId),V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly

def D2(theta, geo, poly=False, forces=False):
    
    if forces==False and poly==False:
        return D1(theta,geo)
    else:
        raise AttributeError("Code not update")
        h=geo.h
        rb=geo.rb
        phi_ie=geo.phi_ie
        phi_e=geo.phi_ie
        phi_o0=geo.phi_o0
        phi_i0=geo.phi_i0
        phi_is=geo.phi_is
        phi_os=geo.phi_os
        ro=rb*(pi-phi_i0+phi_o0)
        
        Nc=getNc(theta,geo=geo)
        (Vd1,dVd1,cxd1,cyd1,fx_pd1,fy_pd1,M_Od1,cd1,V_polyd1,cx_polyd1,cy_polyd1,fxp_polyd1,fyp_polyd1,M_Od1_poly)=D1(theta,geo)
        fx_p=-h*rb*(-sin(theta-phi_e)+(theta+phi_i0-phi_e+2*pi*Nc)*cos(theta-phi_e)+sin(phi_os)-(phi_os-phi_i0+pi)*cos(phi_os))
        fy_p=h*rb*((theta+phi_i0-phi_e+2*pi*Nc)*sin(theta-phi_e)+cos(theta-phi_e)-(-phi_os+phi_i0-pi)*sin(phi_os)+cos(phi_os))
        M_O=-(h*rb**2*(theta-phi_os+2*phi_i0-phi_e+2*pi*Nc-pi)*(theta+phi_os-phi_e+2*pi*Nc+pi))/2
        
        (cx,cy)=(-cxd1+ro*cos(phi_ie-pi/2-theta),-cyd1+ro*sin(phi_ie-pi/2-theta))
    
    if poly==True:
        phi=np.linspace(phi_os+pi,phi_ie-theta-2.0*pi*Nc,1000)
        (xo,yo)=coords_inv(phi, geo, theta, "oi")
        nx=np.zeros_like(phi)
        ny=np.zeros_like(phi)
        (nx,ny)=coords_norm(phi,geo,theta,"oi")
        L=len(xo)
        dA=h*np.sqrt(np.power(xo[1:L]-xo[0:L-1],2)+np.power(yo[1:L]-yo[0:L-1],2))
        fxp_poly=np.sum(dA*(nx[1:L]+nx[0:L-1])/2.0)
        fyp_poly=np.sum(dA*(ny[1:L]+ny[0:L-1])/2.0)
    else:
        fxp_poly=None
        fyp_poly=None
    
    return Vd1,dVd1,cx,cy,fx_p,fy_p,fxp_poly,fyp_poly
    
def DD(theta, geo, poly=False, forces=False):
    
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
    
    if forces==False and poly==False:
        return V, dV
    else:
        raise AttributeError
    if forces==True:
        ################## Force Components #########
        #Arc 1
        fx_p =-hs*geo.ra_arc1*(sin(geo.t2_arc1)-sin(geo.t1_arc1))
        fy_p =+hs*geo.ra_arc1*(cos(geo.t2_arc1)-cos(geo.t1_arc1))
        M_O  =-hs*geo.ra_arc1*((sin(geo.t2_arc1)-sin(geo.t1_arc1))*geo.ya_arc1+(cos(geo.t2_arc1)-cos(geo.t1_arc1))*geo.xa_arc1)
        #Arc 2
        fx_p+=+hs*geo.ra_arc2*(sin(geo.t2_arc2)-sin(geo.t1_arc2))
        fy_p+=-hs*geo.ra_arc2*(cos(geo.t2_arc2)-cos(geo.t1_arc2))
        M_O +=+hs*geo.ra_arc2*((sin(geo.t2_arc2)-sin(geo.t1_arc2))*geo.ya_arc2+(cos(geo.t2_arc2)-cos(geo.t1_arc2))*geo.xa_arc2)
        
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
            M_O+=rx*hs*ny*L-ry*hs*nx*L
        
        #Involute portion
        fx_p+=-hs*(-sin(phi_os)+(phi_os-phi_i0+pi)*cos(phi_os)-sin(phi_is)-(phi_i0-phi_is)*cos(phi_is))*rb
        fy_p+=hs*((-phi_os+phi_i0-pi)*sin(phi_os)-cos(phi_os)-(phi_is-phi_i0)*sin(phi_is)-cos(phi_is))*rb
        M_O +=-(hs*(phi_os-phi_is+pi)*(phi_os+phi_is-2*phi_i0+pi)*rb*rb)/2
    
    if poly==True:
        raise AttributeError('Polygon not coded')
        ##########################################################
        ##                    POLYGON                           ##
        ##########################################################
#        t=np.linspace(geo.t1_arc1,geo.t2_arc1,300)
#        (x_farc1,y_farc1)=(
#            geo.xa_arc1+geo.ra_arc1*cos(t),
#            geo.ya_arc1+geo.ra_arc1*sin(t))
#        (x_oarc1,y_oarc1)=(
#           -geo.xa_arc1-geo.ra_arc1*cos(t)+geo.ro*cos(om),
#           -geo.ya_arc1-geo.ra_arc1*sin(t)+geo.ro*sin(om))
#        (nx_oarc1,ny_oarc1)=(-cos(t),-sin(t))
#        
#        t=np.linspace(geo.t1_arc2,geo.t2_arc2,300)
#        (x_farc2,y_farc2)=(
#            geo.xa_arc2+geo.ra_arc2*cos(t),
#            geo.ya_arc2+geo.ra_arc2*sin(t))
#        (x_oarc2,y_oarc2)=(
#           -geo.xa_arc2-geo.ra_arc2*cos(t)+geo.ro*cos(om),
#           -geo.ya_arc2-geo.ra_arc2*sin(t)+geo.ro*sin(om)) 
#        (nx_oarc2,ny_oarc2)=(+cos(t),+sin(t))
#        
#        phi=np.linspace(phi_is,phi_os+pi,300)
#        (x_finv,y_finv)=coords_inv(phi,geo,theta,'fi')
#        (x_oinv,y_oinv)=coords_inv(phi,geo,theta,'oi')
#        (nx_oinv,ny_oinv)=coords_norm(phi,geo,theta,'oi')
#        
#        x=np.r_[x_farc2[::-1],x_farc1,x_finv,x_oarc2[::-1],x_oarc1,x_oinv,x_farc2[-1]]
#        y=np.r_[y_farc2[::-1],y_farc1,y_finv,y_oarc2[::-1],y_oarc1,y_oinv,y_farc2[-1]]
#        (cx_poly,cy_poly)=polycentroid(x,y)
#        V_poly=geo.h*polyarea(x, y)
#        
#        fxp_poly=0
#        fyp_poly=0
#        MO_poly=0
#        #Arc1
#        L=len(nx_oarc1)
#        dA=hs*np.sqrt(np.power(x_oarc1[1:L]-x_oarc1[0:L-1],2)+np.power(y_oarc1[1:L]-y_oarc1[0:L-1],2))
#        dfxp_poly=dA*(nx_oarc1[1:L]+nx_oarc1[0:L-1])/2.0
#        dfyp_poly=dA*(ny_oarc1[1:L]+ny_oarc1[0:L-1])/2.0
#        fxp_poly=np.sum(dfxp_poly)
#        fyp_poly=np.sum(dfyp_poly)
#        rOx=x_oarc1-geo.ro*cos(phi_e-pi/2-theta)
#        rOx=(rOx[1:L]+rOx[0:L-1])/2
#        rOy=y_oarc1-geo.ro*sin(phi_e-pi/2-theta)
#        rOy=(rOy[1:L]+rOy[0:L-1])/2
#        MO_poly=np.sum(rOx*dfyp_poly-rOy*dfxp_poly)
#        #Arc2
#        L=len(nx_oarc2)
#        dA=hs*np.sqrt(np.power(x_oarc2[1:L]-x_oarc2[0:L-1],2)+np.power(y_oarc2[1:L]-y_oarc2[0:L-1],2))
#        dfxp_poly=dA*(nx_oarc2[1:L]+nx_oarc2[0:L-1])/2.0
#        dfyp_poly=dA*(ny_oarc2[1:L]+ny_oarc2[0:L-1])/2.0
#        fxp_poly+=np.sum(dfxp_poly)
#        fyp_poly+=np.sum(dfyp_poly)
#        rOx=x_oarc2-geo.ro*cos(phi_e-pi/2-theta)
#        rOx=(rOx[1:L]+rOx[0:L-1])/2
#        rOy=y_oarc2-geo.ro*sin(phi_e-pi/2-theta)
#        rOy=(rOy[1:L]+rOy[0:L-1])/2
#        MO_poly+=np.sum(rOx*dfyp_poly-rOy*dfxp_poly)
#        #Involute
#        L=len(y_oinv)
#        dA=hs*np.sqrt(np.power(x_oinv[1:L]-x_oinv[0:L-1],2)+np.power(y_oinv[1:L]-y_oinv[0:L-1],2))
#        dfxp_poly=dA*(nx_oinv[1:L]+nx_oinv[0:L-1])/2.0
#        dfyp_poly=dA*(ny_oinv[1:L]+ny_oinv[0:L-1])/2.0
#        fxp_poly+=np.sum(dfxp_poly)
#        fyp_poly+=np.sum(dfyp_poly)
#        rOx=x_oinv-geo.ro*cos(phi_e-pi/2-theta)
#        rOx=(rOx[1:L]+rOx[0:L-1])/2
#        rOy=y_oinv-geo.ro*sin(phi_e-pi/2-theta)
#        rOy=(rOy[1:L]+rOy[0:L-1])/2
#        MO_poly+=np.sum(rOx*dfyp_poly-rOy*dfxp_poly)
#        #Line
#        x1t=-geo.xa_arc1-geo.ra_arc1*cos(geo.t1_arc1)+ro*cos(om)
#        y1t=-geo.ya_arc1-geo.ra_arc1*sin(geo.t1_arc1)+ro*sin(om)
#        x2t=-geo.xa_arc2-geo.ra_arc2*cos(geo.t1_arc2)+ro*cos(om)
#        y2t=-geo.ya_arc2-geo.ra_arc2*sin(geo.t1_arc2)+ro*sin(om)
#        L=np.sqrt((x2t-x1t)**2+(y2t-y1t)**2)
#        if L>1e-12:
#            Lx=(x2t-x1t)/L
#            Ly=(y2t-y1t)/L
#            nx=-1/np.sqrt(1+Lx**2/Ly**2)
#            ny=Lx/Ly/np.sqrt(1+Lx**2/Ly**2)
#            # Make sure you get the cross product with the normal 
#            # pointing towards the scroll, otherwise flip...
#            if Lx*ny-Ly*nx<0:
#                nx*=-1
#                ny*=-1
#            fxp_poly+=hs*nx*L
#            fyp_poly+=hs*ny*L
#            rx=(x1t+x2t)/2-ro*cos(om)
#            ry=(y2t+y2t)/2-ro*sin(om)
#            MO_poly+=rx*hs*ny*L-ry*hs*nx*L
    else:
        (V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly)=(None,None,None,None,None,None)
        
    return V,dV,cx,cy,fx_p,fy_p,M_O,(V_Oa,dV_Oa,V_Ob,dV_Ob),V_poly,cx_poly,cy_poly,fxp_poly,fyp_poly,MO_poly

def DDD(theta, geo, poly=False, forces=False): 
    
    if poly==False and forces==False:
        V_d1,dV_d1=D1(theta,geo)
        V_d2,dV_d2=D2(theta,geo)
        V_dd,dV_dd=DD(theta,geo)
        V_ddd=V_d1+V_d2+V_dd
        dV_ddd=dV_d1+dV_d2+dV_dd
        return V_ddd,dV_ddd
    else:
        raise AttributeError('Not coded yet')  
 
def phi_s_sa(theta,geo):
    
    h=geo.h
    rb=geo.rb 
    phi_ie=geo.phi_ie
    phi_o0=geo.phi_o0
    phi_i0=geo.phi_i0
    ro=rb*(pi-phi_i0+phi_o0)
    
    b=(-phi_o0+phi_ie-pi);
    D=ro/rb*((phi_i0-phi_ie)*sin(theta)-cos(theta)+1)/(phi_ie-phi_i0);
    B=1.0/2.0*(sqrt(b**2-4.0*D)-b);
    return phi_ie-pi+B-phi_o0
    
def phi_d_dd(theta,geo):
    
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

def Area_d_dd(theta, geo):
    x_fis,y_fis=coords_inv(phi_d_dd(theta,geo),geo,theta,"fi")
    x_oos,y_oos=coords_inv(geo.phi_os,geo,theta,"oo")
    return geo.h*((x_fis-x_oos)**2+(y_fis-y_oos)**2)**0.5

def Area_s_sa(theta,geo):
    return geo.ro*geo.h*(1-cos(theta))
    
if __name__=='__main__':
    print 'This is the base file with scroll geometry.  Running this file doesn\'t do anything'