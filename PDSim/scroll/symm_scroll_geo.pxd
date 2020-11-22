
import numpy as np
cimport numpy as np

import cython
cimport cython

from PDSim.scroll.common_scroll_geo cimport geoVals, HTAnglesClass

from libc.math cimport sqrt,sin,cos,tan,atan2,acos,floor,M_PI as pi,pow,atan

cpdef double fxA(double rb, double phi, double phi0)
cpdef double fyA(double rb, double phi, double phi0)

cpdef double theta_d(geoVals geo) except *

cpdef double radial_leakage_area(double theta, geoVals geo, long key1, long key2, int location = *) except *
cdef _radial_leakage_angles(double theta, geoVals geo, long key1, long key2, double *angle_min, double *angle_max)

cpdef int getNc(double theta, geoVals geo, int path = *) except *

cpdef HTAnglesClass HT_angles(double theta, geoVals geo, key)

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

cpdef CVcoords(CVkey, geoVals geo, double theta)

cpdef double phi_s_sa(double theta, geoVals geo)

@cython.locals(iter=cython.int,phi_os=cython.double,phi_o0=cython.double,phi_ie=cython.double,phi_i0=cython.double,change=cython.double,eps=cython.double,f=cython.double,x1=cython.double,x2=cython.double,x3=cython.double,y1=cython.double,y2=cython.double,phi=cython.double,alpha=cython.double)
cpdef double phi_d_dd(double theta, geoVals geo)

cpdef double Area_d_dd(double theta, geoVals geo)
cpdef double Area_s_sa(double theta, geoVals geo)
cpdef double Area_s_s1_offset(double theta, geoVals geo)