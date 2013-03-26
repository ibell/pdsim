
from PDSim.misc.datatypes cimport arraym
from callbacks cimport StepCallback, HeatTransferCallback

cdef class BaseCycleIntegrator:
    cdef StepCallback step_callback
    
    cpdef store_arrays(self)
    cpdef cycle(self, arraym x0, double min_val, double max_val)