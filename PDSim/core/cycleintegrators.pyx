cdef class BaseCycleIntegrator:
     
    def __init__(self):
        pass
     
    cpdef store_arrays(self):
        raise NotImplementedError("BaseCycleIntegrator is not meant to be called directly")
    
    cpdef cycle(self, arraym x0, double min_val, double max_val):
        raise NotImplementedError("BaseCycleIntegrator is not meant to be called directly")
        