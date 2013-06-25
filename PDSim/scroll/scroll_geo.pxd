
import numpy as np
cimport numpy as np

import cython
cimport cython



from libc.math cimport sqrt,sin,cos,tan,atan2,acos,floor,M_PI as pi,pow

cdef class HTAnglesClass(object):
    cdef public double phi_1_i, phi_2_i, phi_1_o, phi_2_o
    
#Take over the geo class with strict typing
cdef class geoVals:
    cdef public double h,phi_i0,phi_is,phi_ie,phi_e,phi_o0,ro,rb,phi_os,phi_oe,t
    cdef public double xa_arc1,ya_arc1,ra_arc1,t1_arc1,t2_arc1
    cdef public double xa_arc2,ya_arc2,ra_arc2,t1_arc2,t2_arc2
    cdef public double b_line, t1_line, t2_line, m_line
    cdef public double x0_wall,y0_wall,r_wall
    cdef public double delta_radial, delta_flank
    cdef public double phi_ie_offset, delta_suction_offset
    cdef public double cx_scroll, cy_scroll, V_scroll
    cdef public double Vremove
    
cpdef double fxA(double rb, double phi, double phi0)
cpdef double fyA(double rb, double phi, double phi0)

cpdef double theta_d(geoVals geo)
cpdef tuple coords_inv(phi_vec,geoVals geo, double theta, flag = *)
cpdef tuple _coords_inv_np(np.ndarray[np.float_t] phi, geoVals geo, double theta, flag = *)
cpdef tuple _coords_inv_d(double phi, geoVals geo, double theta, flag = *)

cpdef double radial_leakage_area(double theta, geoVals geo, long key1, long key2, int location = *) except *
cdef radial_leakage_angles(double theta, geoVals geo, long key1, long key2, double *angle_min, double *angle_max)

cpdef tuple scroll_wrap(geoVals geo)

cpdef int getNc(double theta, geoVals geo)

cpdef HTAnglesClass HT_angles(double theta, geoVals geo, bytes key)

cpdef tuple SA(double theta, geoVals geo, bint poly=*, bint use_offset = *, double Vremove = *)
cpdef dict SA_forces(double theta, geoVals geo, bint poly = *, bint use_offset = *)

cpdef tuple S1(double theta, geoVals geo, bint poly = *, double theta_0_volume = *, bint use_offset = *)
cpdef dict S1_forces(double theta, geoVals geo, bint poly = *, double theta_0_volume =*, bint use_offset = *)

cpdef tuple S2(double theta, geoVals geo, bint poly = *, double theta_0_volume = *)
cpdef dict S2_forces(double theta, geoVals geo, bint poly = *, double theta_0_volume =*)

cpdef tuple C1(double theta, int alpha, geoVals geo, bint poly=*)
cpdef dict C1_forces(double theta, int alpha, geoVals geo, bint poly = *)

cpdef tuple C2(double theta, int alpha, geoVals geo, bint poly=*)
cpdef dict C2_forces(double theta, int alpha, geoVals geo, bint poly = *)

cpdef tuple D1(double theta, geoVals geo, bint poly=*)
cpdef dict D1_forces(double theta, geoVals geo, bint poly = *)

cpdef tuple D2(double theta, geoVals geo, bint poly=*)
cpdef dict D2_forces(double theta, geoVals geo, bint poly = *)

cpdef tuple DD(double theta, geoVals geo, bint poly=*)
cpdef dict DD_forces(double theta, geoVals geo, bint poly=*)

cpdef tuple DDD(double theta, geoVals geo, bint poly=*) 
cpdef dict DDD_forces(double theta, geoVals geo, bint poly=*) 

cpdef double min2(double, double)
cpdef double max2(double, double)
cpdef double phi_s_sa(double theta, geoVals geo)

@cython.locals(iter=cython.int,phi_os=cython.double,phi_o0=cython.double,phi_ie=cython.double,phi_i0=cython.double,change=cython.double,eps=cython.double,f=cython.double,x1=cython.double,x2=cython.double,x3=cython.double,y1=cython.double,y2=cython.double,phi=cython.double,alpha=cython.double)
cpdef double phi_d_dd(double theta, geoVals geo)

cpdef double Area_d_dd(double theta, geoVals geo)
cpdef double Area_s_sa(double theta, geoVals geo)
cpdef double Area_s_s1_offset(double theta, geoVals geo)