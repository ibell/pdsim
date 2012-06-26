
from CoolProp.State cimport State as StateClass
from CoolProp.State import State as StateClass

cimport numpy as np
import numpy as np

cdef class ControlVolume:
    cdef public StateClass State
    cdef public bint exists
    cdef public object V_dV, V_dV_kwargs
    cdef public bytes key,becomes,discharge_becomes
    
cdef class ControlVolumeCollection(dict):
    cdef dict Idict
    
    cpdef updateStates(self,bytes name1, array1, bytes name2,array2,list keys=?)
    
    @cython.locals(i=cython.int)
    cpdef _T(self)
    