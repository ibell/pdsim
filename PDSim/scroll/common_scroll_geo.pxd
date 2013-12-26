
import numpy as np
cimport numpy as np

import cython
cimport cython

from libc.math cimport sqrt,sin,cos,tan,atan2,acos,floor,M_PI as pi,pow

cdef class VdVstruct:
    cdef public double V,dV

cdef class HTAnglesClass:
    cdef public double phi_1_i, phi_2_i, phi_1_o, phi_2_o, phi_i0, phi_o0
    
cdef enum sides:
    UP
    DOWN
    MID
    
cdef enum compressor_CV_indices:
    keyIsa, keyIs1, keyIs2, keyId1, keyId2,
    keyIdd, keyIddd
    keyIc1_1, keyIc2_1, keyIc1_2, keyIc2_2,
    keyIc1_3, keyIc2_3, keyIc1_4, keyIc2_4,
    keyIc1_5, keyIc2_5, keyIc1_6, keyIc2_6,
    keyIc1_7, keyIc2_7, keyIc1_8, keyIc2_8,
    keyIc1_9, keyIc2_9, keyIc1_10,  keyIc2_10
    
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

cpdef tuple scroll_wrap(geoVals geo)

cpdef tuple coords_norm(phi_vec, geoVals geo, double theta, flag = *)
cpdef tuple coords_inv(phi_vec, geoVals geo, double theta, flag = *)
cpdef tuple _coords_inv_np(np.ndarray[np.float_t] phi, geoVals geo, double theta, flag = *)
cpdef tuple _coords_inv_d(double phi, geoVals geo, double theta, flag = *)

cpdef long get_compressor_CV_index(str key) except *
cpdef long get_compression_chamber_index(long path, long alpha)

cpdef coords_inv_dtheta(phi, geoVals geo, double theta, inv = *)

#ctypedef fused double_or_numpy:
#    cython.double
#    np.array

cpdef double Gr(double phi, geoVals geo, double theta, inv)
cpdef double dGr_dphi(double phi, geoVals geo, double theta, inv)
cpdef double dGr_dtheta(double phi, geoVals geo, double theta, inv)

cpdef double involute_heat_transfer(double hc, double hs, double  rb, 
                                  double phi1, double phi2, double phi0, 
                                  double T_scroll, double T_CV, double dT_dphi, 
                                  double phim)