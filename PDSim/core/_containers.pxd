from CoolProp.State cimport State as StateClass
from PDSim.misc._listmath cimport listm
from libcpp cimport bool

from PDSim.misc.datatypes cimport arraym
from PDSim.misc.datatypes import arraym

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
    
    
cdef class CVArrays(object):
    cdef public arraym T,p,h,rho,V,dV,cp,cv,m,v