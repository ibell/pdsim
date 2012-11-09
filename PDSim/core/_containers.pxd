from CoolProp.State cimport State as StateClass
from PDSim.misc._listmath cimport listm
from libcpp cimport bool

cdef class _ControlVolume:
    cdef public StateClass State
    
cdef class TubeCollection(list):
    cdef dict _Nodes
    
    cpdef dict get_Nodes(self)
    cpdef update(self)
    
cdef class _Tube(object):
    cdef public bytes key1,key2
    cdef public int fixed
    cdef public double Q_add,alpha,L,ID,OD,mdot
    cdef public bool exists
    