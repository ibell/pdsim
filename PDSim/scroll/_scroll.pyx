# cython: profile=True
from PDSim.flow import flow_models
cimport PDSim.flow.flow_models as flow_models

from PDSim.scroll import common_scroll_geo, symm_scroll_geo
cimport PDSim.scroll.common_scroll_geo as common_scroll_geo
cimport PDSim.scroll.symm_scroll_geo as symm_scroll_geo

from libc.math cimport M_PI as pi

cdef int TYPE_RADIAL = flow_models.TYPE_RADIAL
cdef int TYPE_FLANK = flow_models.TYPE_FLANK
cdef int TYPE_DISABLED = flow_models.TYPE_DISABLED

cdef class _Scroll(object):
    
    cpdef dict __cdict__(self):
        return dict(theta = self.theta, 
                    geo = self.geo,
                    HTC = self.HTC)
    
    cpdef double SA_S(self, FlowPath FP):
        FP.A = symm_scroll_geo.Area_s_sa(self.theta, self.geo)
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
        
        FP.A = symm_scroll_geo.radial_leakage_area(self.theta,
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
        
    cpdef double FlankLeakage(self, FlowPath FP, int Ncv_check = -1):
        """
        Calculate the flank leakage flow rate
        
        Parameters
        ----------
        FP : FlowPath
        Ncv_check : int,optional
            If ``Ncv_check`` is greater than -1, this flow path will only be evaluated
            when the number of pairs of compression chambers is equal to this value
        """
        cdef bint _evaluate = False
        cdef double t = -1.0 #Default (not-provided) value
        
        if Ncv_check > -1:
            if Ncv_check == symm_scroll_geo.getNc(self.theta, self.geo):
                _evaluate = True
            else:
                _evaluate = False
        else:
            _evaluate = True
        
        if _evaluate:
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
        else:
            return 0.0
        
