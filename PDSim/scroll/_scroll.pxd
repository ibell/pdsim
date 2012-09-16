from PDSim.flow._flow import _FlowPath
from PDSim.flow._flow cimport _FlowPath

from PDSim.misc._listmath import listm
from PDSim.misc._listmath cimport listm

from PDSim.scroll import scroll_geo
from PDSim.scroll cimport scroll_geo

from PDSim.scroll.scroll_geo import geoVals
from PDSim.scroll.scroll_geo cimport geoVals

cdef class _Scroll(object):
    cdef public geoVals geo
    cdef public double theta
    
    cpdef dict __cdict__(self)
    cpdef double SA_S(self, _FlowPath FP)
    cpdef double Discharge(self,_FlowPath FP)
    cpdef double Inlet_sa(self, _FlowPath FP)
    cpdef double FlankLeakage(self, _FlowPath FP)
    cpdef double RadialLeakage(self, _FlowPath FP, dict kwargs)
    
    cpdef double involute_heat_transfer(self, double hc, double hs, double  rb, 
                                  double phi1, double phi2, double phi0, 
                                  double T_scroll, double T_CV, double dT_dphi, 
                                  double phim)