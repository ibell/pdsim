
import numpy as np
cimport numpy as np

import cython
cimport cython

from libc.math cimport sqrt,sin,cos,tan,atan2,acos,floor,M_PI as pi,pow

cdef class VdVstruct:
    cdef public double V,dV

cdef class HTAnglesClass:
    cdef public double phi_1_i, phi_2_i, phi_1_o, phi_2_o, phi_i0, phi_o0
    
cdef class CVInvolute:
    cdef public double phi_max, phi_min, phi_0
    cdef public int involute
    
cdef class CVInvolutes:
    cdef public CVInvolute Inner, Outer

cdef enum involute_index:
    INVOLUTE_FI
    INVOLUTE_FO
    INVOLUTE_OI
    INVOLUTE_OO
    
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
    cpdef double val_if_symmetric(self, double val)

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

cpdef long get_compressor_CV_index(str key) except *
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

cpdef double involute_heat_transfer(double hc, double hs, double  rb, 
                                  double phi1, double phi2, double phi0, 
                                  double T_scroll, double T_CV, double dT_dphi, 
                                  double phim)