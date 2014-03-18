from PDSim.misc.datatypes import *
from PDSim.misc.datatypes cimport *

########## Heat transfer callbacks ##########
cdef class HeatTransferCallback(object):
    cdef object core 
    cpdef arraym call(self, double t)

cdef class WrappedHeatTransferCallback(HeatTransferCallback):
    cdef object func
    cpdef arraym call(self, double t)
    
########## Lumps energy balance callbacks ##########
cdef class LumpsEnergyBalanceCallback(object):
    cdef object core 
    cpdef arraym call(self)

cdef class WrappedLumpsEnergyBalanceCallback(LumpsEnergyBalanceCallback):
    cdef object func
    cpdef arraym call(self)
    
########## Step callbacks ##########
cdef class StepCallback(object):
    cdef readonly object disable_adaptive
    cdef readonly double h
    cdef object core
    cpdef double call(self, double t, double h, int i) except *
    
cdef class WrappedStepCallback(StepCallback):
    cdef object func
    cpdef double call(self, double t, double h, int i) except *
    
########## Callbacks container ##########    
cdef class CallbackContainer(object):
    cdef public StepCallback step_callback
    cdef public HeatTransferCallback heat_transfer_callback
    cdef public LumpsEnergyBalanceCallback lumps_energy_balance_callback
    cdef public object endcycle_callback