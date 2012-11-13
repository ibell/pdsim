from PDSim.flow._flow import FlowPath
from PDSim.flow._flow cimport FlowPath

from PDSim.scroll import scroll_geo
from PDSim.scroll cimport scroll_geo

from PDSim.scroll.scroll_geo import geoVals
from PDSim.scroll.scroll_geo cimport geoVals
    
cdef class _Scroll(object):
    cdef public geoVals geo
    cdef public double theta
    cdef public double HTC
    
    cpdef dict __cdict__(self)
    cpdef double SA_S(self, FlowPath FP)
    cpdef double Discharge(self,FlowPath FP)
    cpdef double Inlet_sa(self, FlowPath FP)
    cpdef double RadialLeakage(self, FlowPath FP)
    cpdef double FlankLeakage(self, FlowPath FP)
    
    cpdef double calcHT(self, double theta, bytes key, double HTC_tune, double dT_dphi, double phim)
    
    cpdef double involute_heat_transfer(self, double hc, double hs, double  rb, 
                                  double phi1, double phi2, double phi0, 
                                  double T_scroll, double T_CV, double dT_dphi, 
                                  double phim)