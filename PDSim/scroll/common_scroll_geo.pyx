cimport cython
import cython
import numpy as np

cpdef Gr(phi, geoVals geo, double theta, inv):
    
    theta_m = geo.phi_fie - theta + 3.0*pi/2.0
    
    if inv == 'fi':
        return phi*geo.rb**2*(phi**2 - 3*phi*geo.phi_fi0 + 3*geo.phi_fi0**2)/3
    elif inv == 'fo':
        return phi*geo.rb**2*(phi**2 - 3*phi*geo.phi_fo0 + 3*geo.phi_fo0**2)/3
    elif inv == 'oi':
        return geo.rb*(phi**3*geo.rb - 3*phi**2*geo.phi_oi0*geo.rb 
                       + 3*phi*geo.phi_oi0**2*geo.rb 
                       + 3*phi*geo.ro*cos(phi - theta_m) 
                       - 3*geo.phi_oi0*geo.ro*cos(phi - theta_m) 
                       - 3*geo.ro*sin(phi - theta_m))/3
    elif inv == 'oo':
        return geo.rb*(phi**3*geo.rb - 3*phi**2*geo.phi_oo0*geo.rb 
                       + 3*phi*geo.phi_oo0**2*geo.rb 
                       + 3*phi*geo.ro*cos(phi - theta_m) 
                       - 3*geo.phi_oo0*geo.ro*cos(phi - theta_m) 
                       - 3*geo.ro*sin(phi - theta_m))/3
                       
cpdef dGr_dphi(phi, geoVals geo, double theta, inv):
    
    theta_m = geo.phi_fie - theta + 3.0*pi/2.0
    
    if inv == 'fi':
        return geo.rb**2*(phi - geo.phi_fi0)**2
    elif inv == 'fo':
        return geo.rb**2*(phi - geo.phi_fo0)**2
    elif inv == 'oi':
        return geo.rb*(geo.rb*(phi - geo.phi_oi0)**2 - (phi- geo.phi_oi0)*geo.ro*sin(phi - theta_m))
    elif inv == 'oo':
        return geo.rb*(geo.rb*(phi - geo.phi_oo0)**2 - (phi- geo.phi_oo0)*geo.ro*sin(phi - theta_m))

cpdef dGr_dtheta(phi, geoVals geo, double theta, inv):
    
    theta_m = geo.phi_fie - theta + 3.0*pi/2.0
    
    if inv in ['fi','fo']:
        return 0.0
    elif inv == 'oi':
        return geo.rb*geo.ro*(-(phi - geo.phi_oi0)*sin(phi - theta_m) - cos(phi - theta_m))
    elif inv == 'oo':
        return geo.rb*geo.ro*(-(phi - geo.phi_oo0)*sin(phi - theta_m) - cos(phi - theta_m))

cpdef coords_inv_dtheta(phi, geoVals geo, double theta, inv=""):
    """
    Internal function that does the calculation if phi is a double variable 
    """

    rb = geo.rb
    ro = rb*(pi - geo.phi_fi0 + geo.phi_oo0)
    om = geo.phi_fie - theta + 3.0*pi/2.0

    if inv in ["fi", "fo"] :
        dx = 0.0
        dy = 0.0
    elif inv in ["oi", "oo"]:
        dx = +ro*sin(om)
        dy = -ro*cos(om)
    else:
        raise ValueError('flag not valid')
    return (dx, dy)
    
cpdef long get_compressor_CV_index(str key) except *:
    """
    Returns the index defined in the ``compressor_CV_indices`` enum. 
    """
    
    if key == 'sa':
        return keyIsa
    elif key == 's1':
        return keyIs1
    elif key == 's2':
        return keyIs2
    elif key == 'd1':
        return keyId1
    elif key == 'd2':
        return keyId2
    elif key == 'dd':
        return keyIdd
    elif key == 'ddd':
        return keyIddd
    elif key == 'c1.1':
        return keyIc1_1
    elif key == 'c2.1':
        return keyIc2_1
    elif key == 'c1.2':
        return keyIc1_2
    elif key == 'c2.2':
        return keyIc2_2
    elif key == 'c1.3':
        return keyIc1_3
    elif key == 'c2.3':
        return keyIc2_3
    elif key == 'c1.4':
        return keyIc1_4
    elif key == 'c2.4':
        return keyIc2_4
    elif key == 'c1.5':
        return keyIc1_5
    elif key == 'c2.5':
        return keyIc2_5
    elif key == 'c1.6':
        return keyIc1_6
    elif key == 'c2.6':
        return keyIc2_6
    elif key == 'c1.7':
        return keyIc1_7
    elif key == 'c2.7':
        return keyIc2_7
    elif key == 'c1.8':
        return keyIc1_8
    elif key == 'c2.8':
        return keyIc2_8
    elif key == 'c1.9':
        return keyIc1_9
    elif key == 'c2.9':
        return keyIc2_9
    elif key == 'c1.10':
        return keyIc1_10
    elif key == 'c2.10':
        return keyIc2_10
    else:
        return -1
    
cpdef long get_compression_chamber_index(long path, long alpha):
    """
    Return the index for the compression chamber with integers
    
    """
    if path == 1:
        if alpha == 1:
            return keyIc1_1
        elif alpha == 2:
            return keyIc1_2
        elif alpha == 3:
            return keyIc1_3
        elif alpha == 4:
            return keyIc1_4
        elif alpha == 5:
            return keyIc1_5
        elif alpha == 6:
            return keyIc1_6
        elif alpha == 7:
            return keyIc1_7
        elif alpha == 8:
            return keyIc1_8
        elif alpha == 9:
            return keyIc1_9
        elif alpha == 10:
            return keyIc1_10
    elif path == 2:
        if alpha == 1:
            return keyIc2_1
        elif alpha == 2:
            return keyIc2_2
        elif alpha == 3:
            return keyIc2_3
        elif alpha == 4:
            return keyIc2_4
        elif alpha == 5:
            return keyIc2_5
        elif alpha == 6:
            return keyIc2_6
        elif alpha == 7:
            return keyIc2_7
        elif alpha == 8:
            return keyIc2_8
        elif alpha == 9:
            return keyIc2_9
        elif alpha == 10:
            return keyIc2_10
    
# A container for the values for the heat transfer angles
cdef class HTAnglesClass(object):
    def __init__(self):
        pass

#This is a list of all the members in geoVals
geoValsvarlist=['h','ro','rb','t',
                'phi_fi0','phi_fis','phi_fie',
                'phi_fo0','phi_fos','phi_foe',
                'phi_oi0','phi_ois','phi_oie',
                'phi_oo0','phi_oos','phi_ooe',
                'xa_arc1','ya_arc1','ra_arc1','t1_arc1','t2_arc1',
                'xa_arc2','ya_arc2','ra_arc2','t1_arc2','t2_arc2',
                'b_line', 't1_line', 't2_line', 'm_line',
                'x0_wall','y0_wall','r_wall',
                'delta_radial', 'delta_flank',
                'phi_ie_offset','delta_suction_offset',
                'cx_scroll','cy_scroll','V_scroll','Vremove']
 
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
        
    cpdef bint is_symmetric(self):
        """
        Returns true if all the angles for the fixed scroll are the same as for the orbiting scroll
        """
        return (abs(self.phi_fi0-self.phi_oi0) < 1e-14
                and abs(self.phi_fis-self.phi_ois) < 1e-14
                and abs(self.phi_fie-self.phi_oie) < 1e-14
                and abs(self.phi_fo0-self.phi_oo0) < 1e-14
                and abs(self.phi_fos-self.phi_oos) < 1e-14
                and abs(self.phi_foe-self.phi_ooe) < 1e-14
                )
                
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
    rb = geo.rb
    ro = rb*(pi - geo.phi_fi0 + geo.phi_oo0)
    om = geo.phi_fie - theta + 3.0*pi/2.0

    if flag=="fi":
        x = rb*np.cos(phi)+rb*(phi-geo.phi_fi0)*np.sin(phi)
        y = rb*np.sin(phi)-rb*(phi-geo.phi_fi0)*np.cos(phi)
    elif flag=="fo":
        x = rb*np.cos(phi)+rb*(phi-geo.phi_fo0)*np.sin(phi)
        y = rb*np.sin(phi)-rb*(phi-geo.phi_fo0)*np.cos(phi)
    elif flag=="oi":
        x = -rb*np.cos(phi)-rb*(phi-geo.phi_oi0)*np.sin(phi)+ro*np.cos(om)
        y = -rb*np.sin(phi)+rb*(phi-geo.phi_oi0)*np.cos(phi)+ro*np.sin(om)
    elif flag=="oo":
        x = -rb*np.cos(phi)-rb*(phi-geo.phi_oo0)*np.sin(phi)+ro*np.cos(om)
        y = -rb*np.sin(phi)+rb*(phi-geo.phi_oo0)*np.cos(phi)+ro*np.sin(om)
    else:
        raise ValueError('flag not valid')
    return (x,y)
    
cpdef tuple _coords_inv_d(double phi, geoVals geo,double theta, flag=""):
    """
    Internal function that does the calculation if phi is a double variable 
    """

    rb = geo.rb
    ro = rb*(pi - geo.phi_fi0 + geo.phi_oo0)
    om = geo.phi_fie - theta + 3.0*pi/2.0

    if flag=="fi":
        x = rb*cos(phi)+rb*(phi-geo.phi_fi0)*sin(phi)
        y = rb*sin(phi)-rb*(phi-geo.phi_fi0)*cos(phi)
    elif flag=="fo":
        x = rb*cos(phi)+rb*(phi-geo.phi_fo0)*sin(phi)
        y = rb*sin(phi)-rb*(phi-geo.phi_fo0)*cos(phi)
    elif flag=="oi":
        x = -rb*cos(phi)-rb*(phi-geo.phi_oi0)*sin(phi)+ro*cos(om)
        y = -rb*sin(phi)+rb*(phi-geo.phi_oi0)*cos(phi)+ro*sin(om)
    elif flag=="oo":
        x = -rb*cos(phi)-rb*(phi-geo.phi_oo0)*sin(phi)+ro*cos(om)
        y = -rb*sin(phi)+rb*(phi-geo.phi_oo0)*cos(phi)+ro*sin(om)
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

cpdef tuple coords_norm(phi_vec, geoVals geo, double theta,flag="fi"):
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
    
    
cdef double x_antideriv(double phi, double phi_0):
    """ Antiderivative of function for x for scroll_wrap """
    return -cos(phi)*((phi_0-phi)**2-3)-3*sin(phi)*(phi_0-phi)

cdef double y_antideriv(double phi, double phi_0):
    """ Antiderivative of function for y for scroll_wrap """
    return +3*cos(phi)*(phi_0-phi)-sin(phi)*((phi_0-phi)**2-3)
    
@cython.cdivision(True)
cpdef tuple scroll_wrap(geoVals geo):
    """
    Calculate the scroll wrap centroid and volume
    """
        
    # Initial angle of the centerline of the orbiting scroll wrap
    cdef double phi_0 = (geo.phi_oi0 + geo.phi_oo0)/2
    # Ending angle of the centerline of the orbiting scroll wrap
    cdef double phi_e = (geo.phi_oie + geo.phi_ooe)/2
    
    Vwrap = geo.h*geo.rb*geo.t*(phi_e**2/2.0-phi_0**2/2.0-(phi_e-phi_0)*phi_0)
    xstar = geo.h*geo.rb**2*geo.t/Vwrap*(x_antideriv(phi_e, phi_0)-x_antideriv(phi_0, phi_0))
    ystar = geo.h*geo.rb**2*geo.t/Vwrap*(y_antideriv(phi_e, phi_0)-y_antideriv(phi_0, phi_0))
    
    return Vwrap, xstar, ystar