
import numpy as np
cimport numpy as np

import cython
cimport cython

from libc.math cimport sqrt,sin,cos,tan,atan2,acos,floor,M_PI as pi,pow

cdef class VdVstruct:
    """
    A struct with values for volume and derivative of volume w.r.t. crank angle
    """

    cdef public double V
    """ Volume [m^3] """

    cdef public double dV
    """ Derivative of volume with respect to crank angle [m^3/radian] """

cdef class HTAnglesClass:
    """
    A struct with angles associated with the calculations needed for the assesment
    of heat transfer in the chambers of the scroll compressor
    """

    cdef public double phi_1_i
    """ Maximum involute angle on the inner involute of the wrap that forms the outer wall of the CV"""
    
    cdef public double phi_2_i
    """ Minimum involute angle on the inner involute of the wrap that forms the outer wall of the CV"""
        
    cdef public double phi_1_o
    """ Maximum involute angle on the outer involute of the wrap that forms the inner wall of the CV"""
        
    cdef public double phi_2_o
    """ Minimum involute angle on the outer involute of the wrap that forms the inner wall of the CV"""

    cdef public double phi_i0
    """ The initial angle on the inner involute of the wrap that forms the outer wall of the CV"""

    cdef public double phi_o0
    """ The initial angle on the outer involute of the wrap that forms the inner wall of the CV"""

cdef enum involute_index:
    INVOLUTE_FI
    INVOLUTE_FO
    INVOLUTE_OI
    INVOLUTE_OO    

cdef class CVInvolute:
    """
    A simple struct to contain the involute angles
    """

    cdef public double phi_max
    """ The maximum involute angle along this involute """

    cdef public double phi_min
    """ The minimum involute angle along this involute """
    
    cdef public double phi_0
    """ The initial involute angle along this involute """

    cdef public double dphi_max_dtheta
    """ The derivative of phi_max w.r.t. crank angle """
    
    cdef public double dphi_min_dtheta
    """ The derivative of phi_min w.r.t. crank angle """

    cdef public involute_index involute
    """ The involute_index of this involute """
    
cdef class CVInvolutes:
    cdef public CVInvolute Inner
    """ The values for the inner involute of this chamber """
    cdef public CVInvolute Outer
    """ The values for the outer involute of this chamber """
    cdef public bint has_line_1 
    """ Boolean for existence of the line #1 """
    cdef public bint has_line_2
    """ Boolean for existence of the line #2 """
    
cdef enum sides:
    UP
    DOWN
    MID
    
cdef enum compressor_CV_indices:
    keyIsa, keyIs1, keyIs2, keyId1, keyId2,
    keyIdd, keyIddd
    keyIc1_1 = 1001, keyIc2_1 = 2001, 
    keyIc1_2 = 1002, keyIc2_2 = 2002,
    keyIc1_3 = 1003, keyIc2_3 = 2003, 
    keyIc1_4 = 1004, keyIc2_4 = 2004,
    keyIc1_5 = 1005, keyIc2_5 = 2005, 
    keyIc1_6 = 1006, keyIc2_6 = 2006,
    keyIc1_7 = 1007, keyIc2_7 = 2007, 
    keyIc1_8 = 1008, keyIc2_8 = 2008,
    keyIc1_9 = 1009, keyIc2_9 = 2009, 
    keyIc1_10 = 1010,  keyIc2_10 = 2010
    
#Take over the geo class with strict typing
cdef class geoVals:
    cdef public double h,ro,rb,t
    cdef public double phi_fi0,phi_fis,phi_fie,phi_fo0,phi_fos,phi_foe, phi_oi0,phi_ois,phi_oie, phi_oo0,phi_oos,phi_ooe
    cdef public double xa_arc1,ya_arc1,ra_arc1,t1_arc1,t2_arc1
    cdef public double xa_arc2,ya_arc2,ra_arc2,t1_arc2,t2_arc2
    cdef public double b_line, t1_line, t2_line, m_line
    cdef public double x0_wall,y0_wall,r_wall
    cdef public double delta_radial, delta_flank
    cdef public double phi_ie_offset, delta_suction_offset
    cdef public double cx_scroll, cy_scroll, V_scroll
    cdef public double Vremove
    cdef public np.ndarray xvec_disc_port, yvec_disc_port
    
    cpdef bint is_symmetric(self)
    cpdef double val_if_symmetric(self, double val) except *

cpdef bytes involute_index_to_key(int index)

cpdef tuple scroll_wrap(geoVals geo)

ctypedef fused double_or_numpy:
    cython.double
    np.ndarray
    
cpdef tuple coords_norm(phi_vec, geoVals geo, double theta, flag = *)
cpdef tuple coords_inv(phi_vec, geoVals geo, double theta, flag = *)
cpdef tuple _coords_inv_np(np.ndarray[np.float_t] phi, geoVals geo, double theta, flag = *)
cpdef tuple _coords_inv_d(double phi, geoVals geo, double theta, flag = *)

cdef _coords_inv_d_int(double phi, geoVals geo, double theta, int flag, double *x, double *y)
cdef coords_inv_dtheta(double phi, geoVals geo, double theta, int inv, double *dx, double *dy)
cdef _dcoords_inv_dphi_int(double phi, geoVals geo,double theta, int flag, double *dxdphi, double *dydphi)

cpdef long get_compressor_CV_index(object key) except *
cpdef long get_compression_chamber_index(long path, long alpha)

cdef inline bint matchpair(long key1, long key2, long target1, long target2):
    return (key1 == target1 and key2 == target2) or (key2 == target1 and key1 == target2)
    
cpdef inline double min2(double a, double b):
    return a if a<b else b

cpdef inline double max2(double a, double b):
    return a if a>b else b
    
cpdef double_or_numpy plus_one(double_or_numpy x)

cpdef double Gr(double phi, geoVals geo, double theta, int inv)
cpdef double dGr_dphi(double phi, geoVals geo, double theta, int inv)
cpdef double dGr_dtheta(double phi, geoVals geo, double theta, int inv)

cpdef VdVstruct VdV(double theta, geoVals geo, CVInvolutes inv)
cpdef dict forces(double theta, geoVals geo, CVInvolutes inv, double A)
cpdef double fxA(double phi, geoVals geo, double theta, involute_index inv)
cpdef double fyA(double phi, geoVals geo, double theta, involute_index inv)
cpdef double fFx_p(double phi, geoVals geo, double theta, involute_index inv)
cpdef double fFy_p(double phi, geoVals geo, double theta, involute_index inv)
cpdef double fMO_p(double phi, geoVals geo, double theta, involute_index inv)

cpdef double involute_heat_transfer(double hc, double hs, double  rb, 
                                  double phi1, double phi2, double phi0, 
                                  double T_scroll, double T_CV, double dT_dphi, 
                                  double phim)
