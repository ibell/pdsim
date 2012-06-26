from CoolProp.State cimport State as StateClass
from PDSim.misc._listmath cimport listm

cdef class _ControlVolume:
    cdef public StateClass State