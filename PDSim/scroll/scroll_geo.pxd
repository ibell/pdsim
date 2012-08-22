
import cython
cimport cython

import numpy as np
cimport numpy as np

from libc.math cimport sqrt,sin,cos,tan,atan2,acos,floor,M_PI as pi,pow

#Take over the geo class with strict typing
cdef class geoVals:
    cdef public double h,phi_i0,phi_is,phi_ie,phi_e,phi_o0,ro,rb,phi_os,phi_oe,t
    cdef public double xa_arc1,ya_arc1,ra_arc1,t1_arc1,t2_arc1
    cdef public double xa_arc2,ya_arc2,ra_arc2,t1_arc2,t2_arc2
    cdef public double b_line, t1_line, t2_line, m_line
    cdef public double x0_wall,y0_wall,r_wall
    cdef public double delta_radial, delta_flank
    
cpdef double fxA(double rb, double phi, double phi0)
cpdef double fyA(double rb, double phi, double phi0)

cpdef tuple coords_inv(phi_vec,geoVals geo, double theta, flag = *)
cpdef tuple _coords_inv_np(np.ndarray phi, geoVals geo, double theta, flag = *)
cpdef tuple _coords_inv_d(double phi, geoVals geo, double theta, flag = *)

cpdef tuple radial_leakage_angles(double theta, geoVals geo, bytes key1, bytes key2)

cpdef int getNc(double theta, geoVals geo)

cpdef HT_angles(double theta, geoVals geo, bytes key)

cpdef tuple SA(double theta, geoVals geo, bint poly=*, bint forces=*)

cpdef tuple S1(double theta, geoVals geo, bint poly = *, double theta_0_volume = *)
cpdef dict S1_forces(double theta, geoVals geo, bint poly = *, double theta_0_volume =*)

cpdef tuple S2(double theta, geoVals geo, bint poly = *, double theta_0_volume = *)
cpdef dict S2_forces(double theta, geoVals geo, bint poly = *, double theta_0_volume =*)

cpdef tuple C1(double theta, int alpha, geoVals geo, bint poly=*)
cpdef dict C1_forces(double theta, int alpha, geoVals geo, bint poly = *)

cpdef tuple C2(double theta, int alpha, geoVals geo, bint poly=*)
cpdef dict C2_forces(double theta, int alpha, geoVals geo, bint poly = *)

cpdef tuple D1(double theta, geoVals geo, bint poly=*)
cpdef dict D1_forces(double theta, geoVals geo, bint poly = *)

@cython.locals(hs=cython.double,
				xa1=cython.double,
				ya1=cython.double,
				ra1=cython.double,
				ta1_1=cython.double,
				ta1_2=cython.double,
				xa2=cython.double,
				ya2=cython.double,
				ra2=cython.double,
				ta2_1=cython.double,
				ta2_2=cython.double,
				ro=cython.double,
				phi_os=cython.double, 
           		phi_o0=cython.double,
           		phi_is=cython.double, 
           		phi_i0=cython.double,
           		phi_e=cython.double,
           		rb=cython.double,
       			om=cython.double,
           		xoos=cython.double,
   		     	yoos=cython.double,
   		     	V_Oa=cython.double,
                dV_Oa=cython.double,
                x1l=cython.double,
                y1l=cython.double,
                V_Ob=cython.double,
                dV_Ob=cython.double,
                V_Oc=cython.double,
                dV_Oc=cython.double,
                V_Ia=cython.double,
                dV_Ia=cython.double,
                V_Ib=cython.double,
                dV_Ib=cython.double,
                b=cython.double,
                D=cython.double,
                B=cython.double,
                B_prime=cython.double
                )
cpdef DD(double theta, geoVals geo, bint poly=?, bint forces=?)

@cython.locals(V_d1=cython.double,
                dV_d1=cython.double,
                V_d2=cython.double,
                dV_d2=cython.double,
                V_dd=cython.double,
                dV_dd=cython.double,
                V_ddd=cython.double,
                dV_ddd=cython.double)
cpdef DDD(double theta, geoVals geo, bint poly=?, bint forces=?) 

cpdef phi_s_sa(double theta,geoVals geo)

@cython.locals(iter=cython.int,phi_os=cython.double,phi_o0=cython.double,phi_ie=cython.double,phi_i0=cython.double,change=cython.double,eps=cython.double,f=cython.double,x1=cython.double,x2=cython.double,x3=cython.double,y1=cython.double,y2=cython.double,phi=cython.double,alpha=cython.double)
cdef double phi_d_dd(double theta, geoVals geo)

cpdef double Area_d_dd(double theta, geoVals geo)
cpdef double Area_s_sa(double theta, geoVals geo)