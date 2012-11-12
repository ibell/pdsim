
#def rebuildTubeCollection(d):
#    """
#    Used with cPickle to recreate the TubeCollection class
#    """
#    TC=TubeCollection()
#    TC.extend(d)
#    TC.update()
#    print TC.get_Nodes,d,TC,TC[0].TubeFcn
#    return TC


cimport cython

from CoolProp.State import State as StateClass
from CoolProp.State cimport State as StateClass

cdef class TubeCollection(list):
    
    def __init__(self):
        self._Nodes = {}
    
    cpdef dict get_Nodes(self):
        self.update()
        return self._Nodes
    
    property Nodes:
        def __get__(self):
            self.update()
            return self._Nodes

    cpdef update(self):
        """
        _Nodes is a dictionary of flow states for any tubes that exist
        """
        list1=[(Tube.key1,Tube.State1) for Tube in self if Tube.exists]
        list2=[(Tube.key2,Tube.State2) for Tube in self if Tube.exists]
        self._Nodes = dict(list1 + list2)
    
    def __getitem__(self, key):
        for Tube in self:
            if Tube.key1 == key or Tube.key2 == key:
                return Tube
        raise KeyError
    
cdef class CVArrays(object):
    """
    A stub class that contains the arraym arrays of the state variables for
    all the control volumes that are passed into the instantiator
    """
    
    @cython.cdivision(True)
    def __init__(self, CVs, double theta):
        cdef StateClass State
        cdef tuple V_dV
        cdef int N = len(CVs)
        cdef int iCV
        
        self.T = arraym()
        self.T.set_size(N)
        
        self.p = arraym()
        self.p.set_size(N)
        
        self.h = arraym()
        self.h.set_size(N)
        
        self.rho = arraym()
        self.rho.set_size(N)
        
        self.cp = arraym()
        self.cp.set_size(N)
        
        self.cv = arraym()
        self.cv.set_size(N)
        
        self.V = arraym()
        self.V.set_size(N)
        
        self.dV = arraym()
        self.dV.set_size(N)
        
        self.m = arraym()
        self.m.set_size(N)
        
        self.v = arraym()
        self.v.set_size(N)
        
        for iCV in range(N):
            CV = CVs[iCV]
            State = CV.State
            self.T.data[iCV] = State.get_T()
            self.p.data[iCV] = State.get_p()
            self.rho.data[iCV] = State.get_rho()
            self.h.data[iCV] = State.get_h()
            self.cp.data[iCV] = State.get_cp()
            self.cv.data[iCV] = State.get_cv()
            self.V.data[iCV], self.dV.data[iCV] = CV.V_dV(theta,**CV.V_dV_kwargs)
            self.m.data[iCV] = self.rho.data[iCV] * self.V.data[iCV]
            self.v.data[iCV] = 1/self.rho.data[iCV]
    
cdef class _ControlVolume:
    pass

cdef class _Tube:
    pass

cpdef list collect_State_h(list CVList):
    cdef _ControlVolume CV
    return [CV.State.get_h() for CV in CVList]