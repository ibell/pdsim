cdef class LumpsEnergyBalanceCallback(object):
    """
    A wrapper to contain the callback that is called to evaluate the lump energy
    balance, as well as the error on the discharge temperature
    
    Two methods are provided, call() and __call__().  The special method
     __call__() allows an instance of LumpsEnergyBalanceCallback to be called 
     directly like:: 
     
         LEBC = LumpsEnergyBalanceCallback(core)
         LEBC()
         
     or you can call the call() method like::
     
         LEBC.call()
         
     Having both methods allows cython functions to stay at the C++ layer since they
     can call the .call() function with defined argument types and not
     need to come back to python layer for speed 
    """

    def __init__(self, Core):
        self.core = Core
        
    cpdef arraym call(self):
        """ This is the base class function, so don't do anything, return all zeros for the heat transfer rate """
        return NotImplementedError('LumpsEnergyBalanceCallback is not meant to be instantiated directly - you should derive from it')
        
    def __call__(self):
        return self.call()
    
cdef class WrappedLumpsEnergyBalanceCallback(LumpsEnergyBalanceCallback):

    def __init__(self, Core, func):
        self.core = Core
        self.func = func
        
    cpdef arraym call(self):
        """ 
        This function returns the values from the wrapped function
        """
        r = self.func()
        if not isinstance(r,arraym):
            r = arraym(r)
        return r
    
cdef class HeatTransferCallback(object):
    """
    A wrapper to contain the callback that is called to evaluate the heat transfer
    rate for all of the control volumes.  It is not meant to be instantiated
    directly, rather inherited
    
    Two methods are provided, call(t) and __call__(t).  The special method
     __call__(t) allows an instance of StepCallback to be called 
     directly like:: 
     
         SC = StepCallback(core)
         SC(t)
         
     or you can call the call() method like::
     
         SC.call(t)
         
     Having both methods allows cython functions to stay at the C++ layer since they
     can call the .call() function with defined argument types and not
     need to come back to python layer for speed 
    """

    def __init__(self, Core):
        #: The flag which determines whether to enable or disable the adaptive method
        #: If true and the adaptive method is being used, the adaptive solver will
        #: be turned off 
        self.disable_adaptive = False
        self.core = Core
        
    cpdef arraym call(self, double t):
        """ This is the base class function, so don't do anything, return all zeros for the heat transfer rate """
        return NotImplementedError('HeatTransferCallback is not meant to be instantiated directly - you should derive from it')
        
    def __call__(self, double t):
        return self.call(t)
    
cdef class WrappedHeatTransferCallback(HeatTransferCallback):

    def __init__(self, Core, func):
        self.core = Core
        self.func = func
        
    cpdef arraym call(self, double t):
        """ 
        This function returns the values from the wrapped function
        """
        cdef arraym Q = self.func(t)
        if not len(Q) == self.core.CVs.Nexist:
            raise ValueError('Length of Q from wrapped heat transfer callback is not equal to length of number of CV in existence')
        return Q

cdef class StepCallback(object):
    """
    A wrapper to contain the callback that is called when at the beginning of 
    the step, before the step is actually evaluated.
    
    Two methods are provided, call(t,h,i) and __call__(t,h,i).  The special method
     __call__(t,h,i) allows an instance of StepCallback to be called 
     directly like:: 
     
         SC = StepCallback(core)
         SC(t,h,i)
         
     or you can call the call() method like::
     
         SC.call(t,h,i)
         
     Having both methods allows cython functions to stay at the C++ layer since they
     can call the .call() function with defined argument types and not
     need to come back to python layer for speed
    """

    def __init__(self, Core):
        """
        Parameters
        ----------
        Core: The :class:`PDSim.core.core.PDSimCore` subclass
        """
        #: The flag which determines whether to enable or disable the adaptive method
        #: If true and the adaptive method is being used, the adaptive solver will
        #: be turned off 
        self.disable_adaptive = False
        self.core = Core
        
    cpdef double call(self, double t, double h, int i) except *:
        """ This is the base class function, so don't do anything, use the same step size again """
        raise NotImplementedError('StepCallback is not meant to be instantiated directly - you should derive from it')
        
    def __call__(self, double t, double h, int i):
        """ This is the base class function, so don't do anything, use the same step size again """
        return self.call(t,h,i)
    
cdef class WrappedStepCallback(StepCallback):
    """
    This class is intended to provide a python-friendly wrapper of the Cython base class so 
    that high-level python code can seamlessly interface through into Cython.  It is not as fast
    as developing a cython Callback function, but it is easier to do
    """
    
    def __init__(self, Core, func):
        """
        Parameters
        ----------
        Core: The :class:`PDSim.core.core.PDSimCore` subclass
        func: The function that will be called
        
        """
        StepCallback.__init__(self,Core)
        self.func = func
        
    cpdef double call(self, double t, double h, int i) except *:
        """ 
        This function returns the values from the wrapped function
        """
        vals = self.func(t,h,i)
        try:
            self.disable_adaptive,h = vals
        except TypeError:
            raise TypeError('step_callback must return a bool,float pair, returned the values:'+str(vals))
        return h
    
cdef class CallbackContainer(object):
    # A cython structure to be used to hold references to all the callbacks
    pass