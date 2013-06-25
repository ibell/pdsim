# cython: profile=True
from PDSim.flow import flow_models
from PDSim.flow cimport flow_models

from PDSim.scroll import scroll_geo
from PDSim.scroll cimport scroll_geo

from PDSim.scroll.scroll_geo import Area_s_sa
from PDSim.scroll.scroll_geo cimport Area_s_sa

from libc.math cimport M_PI as pi

cdef int TYPE_RADIAL = flow_models.TYPE_RADIAL
cdef int TYPE_FLANK = flow_models.TYPE_FLANK

cdef class _Scroll(object):
    
    cpdef dict __cdict__(self):
        return dict(theta = self.theta, 
                    geo = self.geo,
                    HTC = self.HTC)
    
    cpdef double SA_S(self, FlowPath FP):
        FP.A=scroll_geo.Area_s_sa(self.theta,self.geo)
        try:
            return flow_models.IsentropicNozzle(FP.A,FP.State_up,FP.State_down)
        except ZeroDivisionError:
            return 0.0  
        
    cpdef double Discharge(self,FlowPath FP):
        FP.A=pi*0.01**2/4.0
        try:
            return flow_models.IsentropicNozzle(FP.A,FP.State_up,FP.State_down)
        except ZeroDivisionError:
            return 0.0
        
    cpdef double Inlet_sa(self, FlowPath FP):
        FP.A=pi*0.03**2/4.0
        return flow_models.IsentropicNozzle(FP.A,FP.State_up,FP.State_down)
        
    cpdef double RadialLeakage(self, FlowPath FP, double t = -1):
        """
        Calculate the radial leakage flow rate
        
        Parameters
        ----------
        FP : FlowPath
            
        t : float,optional
            The thickness of the wrap to be used.  If not provided, the scroll
            wrap thickness
        """
        #Calculate the area
        #Arc length of the upstream part of the flow path
        
        FP.A = scroll_geo.radial_leakage_area(self.theta,
                                              self.geo,
                                              FP.key1Index,
                                              FP.key2Index)
        
        if FP.A == 0.0:
            return 0.0
        
        # Allow you to change the length for the radial leakage path
        # by passing in a length other than the thickness of the scroll wrap
        # 
        if t <= 0:
            t = self.geo.t
        
        return flow_models.FrictionCorrectedIsentropicNozzle(FP.A, 
                                                             FP.State_up, 
                                                             FP.State_down, 
                                                             self.geo.delta_radial, 
                                                             TYPE_RADIAL, 
                                                             t)
        
    cpdef double FlankLeakage(self,FlowPath FP):
        """
        Calculate the flank leakage flow rate
        """
        cdef double t = -1.0 #Default (not-provided) value
        #Calculate the area
        FP.A = self.geo.h*self.geo.delta_flank
        return flow_models.FrictionCorrectedIsentropicNozzle(
                             FP.A,
                             FP.State_up,
                             FP.State_down,
                             self.geo.delta_flank,
                             TYPE_FLANK,
                             t,
                             self.geo.ro
                             )
        
    cpdef double calcHT(self, double theta, bytes key, double HTC_tune, double dT_dphi, double phim):
        cdef scroll_geo.HTAnglesClass angles
        cdef double hc,Q_outer_wrap,Q_inner_wrap,Q_d1,Q_d2,Q_dd,T_CV,T_scroll
        
        ## If HT is turned off, no heat transfer
        if HTC_tune <= 0.0 or key.startswith('inj') or key == 'sa' or key == 'dd':
            return 0.0
        elif key == 'ddd':
            # ddd is a combination of the heat transfer in the d1, d2, and
            # dd chambers
            Q_d1 = self.calcHT(theta,str('d1'),HTC_tune,dT_dphi,phim)
            Q_d2 = self.calcHT(theta,str('d2'),HTC_tune,dT_dphi,phim)
            Q_d3 = self.calcHT(theta,str('dd'),HTC_tune,dT_dphi,phim)
                
        #print 'calculate HTC'
        #TODO: calculate HTC
        hc = self.HTC #[kW/m2/K]
            
        #Get the bounding angles for the control volume
        angles = scroll_geo.HT_angles(theta, self.geo, key)
        
        T_scroll = self.Tlumps[0]
        T_CV = self.CVs[key].State.T
        # The heat transfer rate of the inner involute on 
        # the outer wrap of the chamber
        Q_outer_wrap = self.involute_heat_transfer(hc, 
                                                   self.geo.h, 
                                                   self.geo.rb, 
                                                   angles.phi_1_i, 
                                                   angles.phi_2_i, 
                                                   self.geo.phi_i0, 
                                                   T_scroll,
                                                   T_CV, 
                                                   dT_dphi, 
                                                   phim)
        
        # The heat transfer rate of the outer involute on 
        # the inner wrap of the chamber
        Q_inner_wrap = self.involute_heat_transfer(hc, 
                                                   self.geo.h, 
                                                   self.geo.rb, 
                                                   angles.phi_1_o, 
                                                   angles.phi_2_o, 
                                                   self.geo.phi_o0,
                                                   T_scroll,
                                                   T_CV, 
                                                   dT_dphi, 
                                                   phim)
        
        return HTC_tune *(Q_outer_wrap + Q_inner_wrap)
        
            
                
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