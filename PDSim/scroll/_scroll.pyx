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
        
    cpdef double involute_heat_transfer(self, double hc, double hs, double  rb, 
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