from PDSim.flow import flow_models
from PDSim.flow cimport flow_models

from PDSim.scroll import scroll_geo
from PDSim.scroll cimport scroll_geo

from PDSim.scroll.scroll_geo import Area_s_sa
from PDSim.scroll.scroll_geo cimport Area_s_sa

from libc.math cimport M_PI as pi

cdef class _Scroll(object):
    
    cpdef dict __cdict__(self):
        return dict(theta = self.theta, geo = self.geo)
    
    cpdef double SA_S(self, _FlowPath FP):
        FP.A=scroll_geo.Area_s_sa(self.theta,self.geo)
        try:
            return flow_models.IsentropicNozzle(FP.A,FP.State_up,FP.State_down)
        except ZeroDivisionError:
            return 0.0  
        
    cpdef double Discharge(self,_FlowPath FP):
        FP.A=pi*0.01**2/4.0
        try:
            return flow_models.IsentropicNozzle(FP.A,FP.State_up,FP.State_down)
        except ZeroDivisionError:
            return 0.0
        
    cpdef double Inlet_sa(self, _FlowPath FP):
        FP.A=pi*0.03**2/4.0
        return flow_models.IsentropicNozzle(FP.A,FP.State_up,FP.State_down)
    
    cpdef double FlankLeakage(self, _FlowPath FP):
        #Calculate the area
        FP.A=self.geo.h*self.geo.delta_flank
        try:
            return flow_models.IsentropicNozzle(FP.A,FP.State_up,FP.State_down)
        except ZeroDivisionError:
            return 0.0