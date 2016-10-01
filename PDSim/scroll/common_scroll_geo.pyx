from __future__ import division
cimport cython
import cython
import numpy as np

cpdef double_or_numpy plus_one(double_or_numpy x):
    return x + 1.0

cdef class VdVstruct:
    def __init__(self, V, dV):
        self.V = V
        self.dV = dV
        
    def __repr__(self):
        return str(list([self.V,self.dV]))

cdef class HTAnglesClass:
    def __repr__(self):
        s = ''
        for k in ['phi_1_i','phi_2_i','phi_1_o','phi_2_o','phi_i0','phi_o0']:
            s += k + ' : ' + str(getattr(self,k)) + '\n'
        return s
        
cdef class CVInvolute:
    def __init__(self):
        pass
        
cdef class CVInvolutes:
    def __init__(self):
        self.Inner = CVInvolute.__new__(CVInvolute)
        self.Outer = CVInvolute.__new__(CVInvolute)
        
    def __repr__(self):
        s = ''
        s += "Outer.involute = {i:s}\n".format(i=involute_index_to_key(self.Outer.involute))
        s += "Outer.phi_0 = {i:g}\n".format(i=self.Outer.phi_0)
        s += "Outer.phi_max = {i:g}\n".format(i=self.Outer.phi_max)
        s += "Outer.phi_min = {i:g}\n".format(i=self.Outer.phi_min)
        s += "Inner.involute = {i:s}\n".format(i=involute_index_to_key(self.Inner.involute))
        s += "Inner.phi_0 = {i:g}\n".format(i=self.Inner.phi_0)
        s += "Inner.phi_max = {i:g}\n".format(i=self.Inner.phi_max)
        s += "Inner.phi_min = {i:g}\n".format(i=self.Inner.phi_min)
        s += "has_line_1 = {i:g}\n".format(i=self.has_line_1)
        s += "has_line_2 = {i:g}".format(i=self.has_line_2)
        return s

cpdef VdV_common(double theta, geoVals geo, CVInvolutes inv):
    """
    Evaluate V and dV/dtheta in a generalized manner for a chamber
    """
    cdef double A_i, A_o, A_line_1 = 0, A_line_2 = 0, x_1, x_2, y_1, y_2
    cdef double dA_line_1_dtheta = 0, dA_line_2_dtheta = 0, dx_1_dtheta, dy_1_dtheta, dx_2_dtheta, dy_2_dtheta

    ## ------------------------ VOLUME -------------------------------

    A_i = 0.5*(Gr(inv.Outer.phi_max, geo, theta, inv.Outer.involute) - Gr(inv.Outer.phi_min, geo, theta, inv.Outer.involute))
    if inv.has_line_1:
        _coords_inv_d_int(inv.Outer.phi_max, geo, theta, inv.Outer.involute, &x_1, &y_1)
        _coords_inv_d_int(inv.Inner.phi_max, geo, theta, inv.Inner.involute, &x_2, &y_2)
        A_line_1 = 0.5*(x_1*y_2 - x_2*y_1)
    A_o = 0.5*(Gr(inv.Inner.phi_min, geo, theta, inv.Inner.involute) - Gr(inv.Inner.phi_max, geo, theta, inv.Inner.involute))
    if inv.has_line_2:
        _coords_inv_d_int(inv.Inner.phi_min, geo, theta, inv.Inner.involute, &x_1, &y_1)
        _coords_inv_d_int(inv.Outer.phi_min, geo, theta, inv.Outer.involute, &x_2, &y_2)
        A_line_2 = 0.5*(x_1*y_2 - x_2*y_1)
    
    V = geo.h*(A_i + A_line_1 + A_o + A_line_2)
    
    ## ------------------------ DERIVATIVE -------------------------------
    
    dA_i_dtheta = 0.5*(dGr_dphi(inv.Outer.phi_max, geo, theta, inv.Outer.involute)*inv.Outer.dphi_max_dtheta
                      +dGr_dtheta(inv.Outer.phi_max, geo, theta, inv.Outer.involute)
                      -dGr_dphi(inv.Outer.phi_min, geo, theta, inv.Outer.involute)*inv.Outer.dphi_min_dtheta
                      -dGr_dtheta(inv.Outer.phi_min, geo, theta, inv.Outer.involute)
                      )

    if inv.has_line_1:
        coords_inv_dtheta(inv.Outer.phi_max, geo, theta, inv.Outer.involute, &dx_1_dtheta, &dy_1_dtheta)
        coords_inv_dtheta(inv.Inner.phi_max, geo, theta, inv.Inner.involute, &dx_2_dtheta, &dy_2_dtheta)
        dA_line_1_dtheta = 0.5*(x_1*dy_2_dtheta + y_2*dx_1_dtheta - x_2*dy_1_dtheta - y_1*dx_2_dtheta)
    
    dA_o_dtheta = 0.5*(dGr_dphi(inv.Inner.phi_min, geo, theta, inv.Inner.involute)*inv.Inner.dphi_min_dtheta
                       +dGr_dtheta(inv.Inner.phi_min, geo, theta, inv.Inner.involute)
                       -dGr_dphi(inv.Inner.phi_max, geo, theta, inv.Inner.involute)*inv.Inner.dphi_max_dtheta
                       -dGr_dtheta(inv.Inner.phi_max, geo, theta, inv.Inner.involute)
                       )

    if inv.has_line_2:
        coords_inv_dtheta(inv.Inner.phi_min, geo, theta, inv.Inner.involute, &dx_1_dtheta, &dy_1_dtheta)
        coords_inv_dtheta(inv.Outer.phi_min, geo, theta, inv.Outer.involute, &dx_2_dtheta, &dy_2_dtheta)
        dA_line_2_dtheta = 0.5*(x_1*dy_2_dtheta + y_2*dx_1_dtheta - x_2*dy_1_dtheta - y_1*dx_2_dtheta)

    dV = geo.h*(dA_i_dtheta + dA_line_1_dtheta + dA_o_dtheta + dA_line_2_dtheta)
    
    cdef VdVstruct VdV = VdVstruct.__new__(VdVstruct)
    VdV.V = V
    VdV.dV = dV
    return VdV

cpdef double fxA(double phi, geoVals geo, double theta, involute_index inv):
    """
    The anti-derivative for the x-coordinate of the centroid
    
    Parameters
    ----------
    phi : float
        The involute angle
    geo : geoVals
        The structure with the geometry obtained from get_geo()
    theta : float
        The crank angle, between 0 and 2*pi
    inv : involute_index
        The key for the involute to be considered
    """
    
    cdef double Theta = geo.phi_fie - theta - pi/2, r_b = geo.rb, r_o = geo.ro
    
    if inv == INVOLUTE_FI:
        phi_0 = geo.phi_fi0
        return r_b**3*(-3*(phi - phi_0)*((phi - phi_0)**2 - 6)*cos(phi) + (phi - phi_0)*((phi - phi_0)**2 - 3)*cos(phi)**3 
                      + 3*(2*(phi - phi_0)**2 - 5)*sin(phi) + (3*(phi - phi_0)**2 - 1)*sin(phi)**3)/6
    elif inv == INVOLUTE_FO:
        phi_0 = geo.phi_fo0
        return r_b**3*(-3*(phi - phi_0)*((phi - phi_0)**2 - 6)*cos(phi) + (phi - phi_0)*((phi - phi_0)**2 - 3)*cos(phi)**3 
                      + 3*(2*(phi - phi_0)**2 - 5)*sin(phi) + (3*(phi - phi_0)**2 - 1)*sin(phi)**3)/6
    elif inv == INVOLUTE_OO:
        phi_0 = geo.phi_oo0 
    elif inv == INVOLUTE_OI:
        phi_0 = geo.phi_oi0
    else:
        raise ValueError("Involute is invalid")
    
    return -r_b*(-phi*r_b*r_o*(phi**2 - 3*phi*phi_0 + 3*phi_0**2 + 3)*cos(Theta)/3 
                + r_b**2*(phi - phi_0)*((phi - phi_0)**2 - 3)*cos(phi)**3/3 
                + r_b**2*(3*(phi - phi_0)**2 - 1)*sin(phi)**3/3 + 2*r_b*r_o*(phi - phi_0)*cos(Theta)*cos(phi)**2 
                - (phi - phi_0)*(r_b**2*((phi - phi_0)**2 - 6) + r_o**2*cos(Theta)**2)*cos(phi) 
                + (r_b**2*(2*(phi - phi_0)**2 - 5) + r_b*r_o*((phi - phi_0)**2 - 1)*cos(Theta)*cos(phi) + r_o**2*cos(Theta)**2)*sin(phi)
                )/2
                        
cpdef double fyA(double phi, geoVals geo, double theta, involute_index inv):
    """
    The anti-derivative for the y-coordinate of the centroid
    
    Parameters
    ----------
    phi : float
        The involute angle
    geo : geoVals
        The structure with the geometry obtained from get_geo()
    theta : float
        The crank angle, between 0 and 2*pi
    inv : involute_index
        The key for the involute to be considered
    """
    
    cdef double Theta = geo.phi_fie - theta - pi/2, r_b = geo.rb, r_o = geo.ro
    
    if inv == INVOLUTE_FI:
        phi_0 = geo.phi_fi0
        return r_b**3*(-(phi - phi_0)*((phi - phi_0)**2 - 3)*sin(phi)**3 + 6*(phi - phi_0)*sin(phi) 
                       - 3*(2*(phi - phi_0)**2 - 3)*cos(phi) + (3*(phi - phi_0)**2 - 1)*cos(phi)**3)/3
    elif inv == INVOLUTE_FO:
        phi_0 = geo.phi_fo0
        return r_b**3*(-(phi - phi_0)*((phi - phi_0)**2 - 3)*sin(phi)**3 + 6*(phi - phi_0)*sin(phi) 
                      - 3*(2*(phi - phi_0)**2 - 3)*cos(phi) + (3*(phi - phi_0)**2 - 1)*cos(phi)**3)/3
    elif inv == INVOLUTE_OO:
        phi_0 = geo.phi_oo0
    elif inv == INVOLUTE_OI:
        phi_0 = geo.phi_oi0
    else:
        raise ValueError("Involute is invalid")

    return -r_b*(-r_b**2*(phi - phi_0)*((phi - phi_0)**2 - 3)*sin(phi)**3/3 - r_b**2*(2*(phi - phi_0)**2 - 3)*cos(phi) 
    + r_b**2*(3*(phi - phi_0)**2 - 1)*cos(phi)**3/3 - r_b*r_o*(phi - phi_0 - 1)*(phi - phi_0 + 1)*cos(Theta)/2 
    + r_b*(4*r_b*(phi - phi_0) + 2*r_o*(phi - phi_0)*cos(Theta + phi) + r_o*((phi - phi_0)**2 - 1)*sin(Theta + phi))*sin(phi)/2 
    - r_o*(r_b*(phi**3 - 3*phi**2*phi_0 + 3*phi*phi_0**2 - 3*phi + 6*phi_0) + 6*r_o*(phi - phi_0)*cos(Theta)*cos(phi) - 6*r_o*sin(phi)*cos(Theta))*sin(Theta)/6)

cpdef double fFx_p(double phi, geoVals geo, double theta, involute_index inv):
    """
    The anti-derivative of the x-coordinate force term (Table 2)
    
    Parameters
    ----------
    phi : float
        The involute angle
    geo : geoVals
        The structure with the geometry obtained from get_geo()
    theta : float
        The crank angle, between 0 and 2*pi
    inv : involute_index
        The key for the involute to be considered
    """
    if inv == INVOLUTE_FI:
        return geo.h*geo.rb*((geo.phi_fi0-phi)*cos(phi)+sin(phi))
    elif inv == INVOLUTE_OO:
        return geo.h*geo.rb*((geo.phi_oo0-phi)*cos(phi)+sin(phi))
    elif inv == INVOLUTE_FO:
        return geo.h*geo.rb*((phi-geo.phi_fo0)*cos(phi)-sin(phi))
    elif inv == INVOLUTE_OI:
        return geo.h*geo.rb*((phi-geo.phi_oi0)*cos(phi)-sin(phi))
    else:
        raise ValueError("Involute is invalid")

cpdef double fFy_p(double phi, geoVals geo, double theta, involute_index inv):
    """
    The anti-derivative of the y-coordinate force term (Table 2)
    
    Parameters
    ----------
    phi : float
        The involute angle
    geo : geoVals
        The structure with the geometry obtained from get_geo()
    theta : float
        The crank angle, between 0 and 2*pi
    inv : involute_index
        The key for the involute to be considered
    """
    if inv == INVOLUTE_FI:
        return geo.h*geo.rb*((geo.phi_fi0-phi)*sin(phi)-cos(phi))
    elif inv == INVOLUTE_OO:
        return geo.h*geo.rb*((geo.phi_oo0-phi)*sin(phi)-cos(phi))
    elif inv == INVOLUTE_FO:
        return geo.h*geo.rb*((phi-geo.phi_fo0)*sin(phi)+cos(phi))
    elif inv == INVOLUTE_OI:
        return geo.h*geo.rb*((phi-geo.phi_oi0)*sin(phi)+cos(phi))
    else:
        raise ValueError("Involute is invalid")
        
cpdef double fMO_p(double phi, geoVals geo, double theta, involute_index inv):
    """
    The anti-derivative of the spinning moment calculation (Table 2)
    
    Parameters
    ----------
    phi : float
        The involute angle
    geo : geoVals
        The structure with the geometry obtained from get_geo()
    theta : float
        The crank angle, between 0 and 2*pi
    inv : involute_index
        The key for the involute to be considered
    """
    if inv == INVOLUTE_OO:
        return geo.h*geo.rb**2*phi/2.0*(phi-2*geo.phi_oo0)
    elif inv == INVOLUTE_OI:
        return -geo.h*geo.rb**2*phi/2.0*(phi-2*geo.phi_oi0)
    else:
        return 1e99
        
cpdef bytes involute_index_to_key(int index):
    """
    Return the string associated with a given index from the common_scroll_geo.involute_index enumeration
    """
    
    if index == INVOLUTE_FI:
        return bytes('fi')
    elif index == INVOLUTE_FO:
        return bytes('fo')
    elif index == INVOLUTE_OI:
        return bytes('oi')
    elif index == INVOLUTE_OO:
        return bytes('oo')
    else:
        return bytes('')
    
cpdef double Gr(double phi, geoVals geo, double theta, int inv):
    """
    The antiderivative of the area integration term, where 
    
    .. math::
    
        Gr \\equiv \\int\\left[\\left(-y\\frac{dx(\\phi)}{d\\phi}+x\\frac{dy(\\phi)}{d\\phi}\\right)d\\phi\\right]
    """
    theta_m = geo.phi_fie - theta + 3.0*pi/2.0
    if inv == INVOLUTE_FI:
        return phi*geo.rb**2*(phi**2 - 3*phi*geo.phi_fi0 + 3*geo.phi_fi0**2)/3
    elif inv == INVOLUTE_FO:
        return phi*geo.rb**2*(phi**2 - 3*phi*geo.phi_fo0 + 3*geo.phi_fo0**2)/3
    elif inv == INVOLUTE_OI:
        return geo.rb*(phi**3*geo.rb - 3*phi**2*geo.phi_oi0*geo.rb 
                       + 3*phi*geo.phi_oi0**2*geo.rb 
                       + 3*(phi-geo.phi_oi0)*geo.ro*cos(phi - theta_m)
                       - 3*geo.ro*sin(phi - theta_m))/3
    elif inv == INVOLUTE_OO:
        return geo.rb*(phi**3*geo.rb - 3*phi**2*geo.phi_oo0*geo.rb 
                       + 3*phi*geo.phi_oo0**2*geo.rb 
                       + 3*(phi-geo.phi_oo0)*geo.ro*cos(phi - theta_m) 
                       - 3*geo.ro*sin(phi - theta_m))/3
                       
cpdef double dGr_dphi(double phi, geoVals geo, double theta, int inv):
    """
    The partial derivative of Gr with respect to phi with theta held constant
    """
    
    THETA = geo.phi_fie - theta - pi/2.0
    
    if inv == INVOLUTE_FI:
        return geo.rb**2*(phi - geo.phi_fi0)**2
    elif inv == INVOLUTE_FO:
        return geo.rb**2*(phi - geo.phi_fo0)**2
    elif inv == INVOLUTE_OI:
        return geo.rb*(geo.rb*(phi - geo.phi_oi0)**2 + (phi- geo.phi_oi0)*geo.ro*sin(THETA - phi))
    elif inv == INVOLUTE_OO:
        return geo.rb*(geo.rb*(phi - geo.phi_oo0)**2 + (phi- geo.phi_oo0)*geo.ro*sin(THETA - phi))

cpdef double dGr_dtheta(double phi, geoVals geo, double theta, int inv):
    """
    The partial derivative of Gr with respect to theta with phi held constant
    """
    
    THETA = geo.phi_fie - theta - pi/2.0
    
    if inv == INVOLUTE_FI or inv == INVOLUTE_FO:
        return 0.0
    elif inv == INVOLUTE_OI:
        return geo.rb*geo.ro*((phi - geo.phi_oi0)*sin(THETA - phi) - cos(THETA - phi))
    elif inv == INVOLUTE_OO:
        return geo.rb*geo.ro*((phi - geo.phi_oo0)*sin(THETA - phi) - cos(THETA - phi))

cdef coords_inv_dtheta(double phi, geoVals geo, double theta, int inv, double *dx, double *dy):
    """
    Internal function that does the calculation if phi is a double variable 
    """

    rb = geo.rb
    ro = rb*(pi - geo.phi_fi0 + geo.phi_oo0)
    om = geo.phi_fie - theta + 3.0*pi/2.0

    if inv == INVOLUTE_FI or inv == INVOLUTE_FO:
        dx[0] = 0.0
        dy[0] = 0.0
    elif inv == INVOLUTE_OI or inv == INVOLUTE_OO:
        dx[0] = +ro*sin(om)
        dy[0] = -ro*cos(om)
    else:
        raise ValueError('flag not valid')
    
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
    return 1000*path+alpha

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
        
    cpdef double val_if_symmetric(self, double val):
        """ Returns the value if symmetric, throws ValueError otherwise """
        if self.is_symmetric():
            return val
        else:
            raise ValueError('Symmetric angle requested (phi_ie, phi_is, phi_i0, etc.) but the geometry is not symmetric')
            
    property phi_ie:
        """ Inner Ending Angle """
        def __get__(self):
            return self.val_if_symmetric(self.phi_fie)        
    
    property phi_is:
        """ Inner Starting Angle """
        def __get__(self):
            return self.val_if_symmetric(self.phi_fis)
            
    property phi_i0:
        """ Inner Initial Angle """
        def __get__(self):
            return self.val_if_symmetric(self.phi_fi0)
            
    property phi_oe:
        """ Outer Ending Angle """
        def __get__(self):
            return self.val_if_symmetric(self.phi_foe)        
    
    property phi_os:
        """ Outer Starting Angle """
        def __get__(self):
            return self.val_if_symmetric(self.phi_fos)
            
    property phi_o0:
        """ Outer Initial Angle """
        def __get__(self):
            return self.val_if_symmetric(self.phi_fo0)        
    
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
    
cdef _coords_inv_d_int(double phi, geoVals geo,double theta, int flag, double *x, double *y):
    """
    Internal function that does the calculation if phi is a double variable 
    """

    rb = geo.rb
    ro = rb*(pi - geo.phi_fi0 + geo.phi_oo0)
    om = geo.phi_fie - theta + 3.0*pi/2.0

    if flag == INVOLUTE_FI:
        x[0] = rb*cos(phi)+rb*(phi-geo.phi_fi0)*sin(phi)
        y[0] = rb*sin(phi)-rb*(phi-geo.phi_fi0)*cos(phi)
    elif flag == INVOLUTE_FO:
        x[0] = rb*cos(phi)+rb*(phi-geo.phi_fo0)*sin(phi)
        y[0] = rb*sin(phi)-rb*(phi-geo.phi_fo0)*cos(phi)
    elif flag == INVOLUTE_OI:
        x[0] = -rb*cos(phi)-rb*(phi-geo.phi_oi0)*sin(phi)+ro*cos(om)
        y[0] = -rb*sin(phi)+rb*(phi-geo.phi_oi0)*cos(phi)+ro*sin(om)
    elif flag == INVOLUTE_OO:
        x[0] = -rb*cos(phi)-rb*(phi-geo.phi_oo0)*sin(phi)+ro*cos(om)
        y[0] = -rb*sin(phi)+rb*(phi-geo.phi_oo0)*cos(phi)+ro*sin(om)
    else:
        raise ValueError('flag not valid')
    
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
    
cpdef double involute_heat_transfer(double hc, double hs, double  rb, 
                                  double phi1, double phi2, double phi0, 
                                  double T_scroll, double T_CV, double dT_dphi, 
                                  double phim):
        """
        This function evaluates the anti-derivative of 
        the differential of involute heat transfer, and returns the amount of scroll-
        wall heat transfer in kW
        
        Parameters
        ----------
        hc : float
            Heat transfer coefficient [kW/m2/K]
        hs : float
            Scroll wrap height [m]
        rb : float
            Base circle radius [m]
        phi1 : float
            Larger involute angle [rad]
        phi2 : float
            Smaller involute angle [rad]
        phi0 : float
            Initial involute angle [rad]
        T_scroll : float
            Lump temperature of the scroll wrap [K]
        T_CV : float
            Temperature of the gas in the CV [K]
        dT_dphi : float
            Derivative of the temperature along the scroll wrap [K/rad]
        phim : float
            Mean involute angle of wrap used for heat transfer [rad]
        
        Notes
        -----
        ``phi1`` and ``phi2`` are defined such that ``phi1`` is always the
        larger involute angle in value
        """
        term1=hc*hs*rb*( (phi1*phi1/2.0-phi0*phi1)*(T_scroll-T_CV)
            +dT_dphi*(phi1*phi1*phi1/3.0-(phi0+phim)*phi1*phi1/2.0+phi0*phim*phi1))
        term2=hc*hs*rb*( (phi2*phi2/2.0-phi0*phi2)*(T_scroll-T_CV)
            +dT_dphi*(phi2*phi2*phi2/3.0-(phi0+phim)*phi2*phi2/2.0+phi0*phim*phi2))
        return term1-term2;