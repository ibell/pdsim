
#def rebuildTubeCollection(d):
#    """
#    Used with cPickle to recreate the TubeCollection class
#    """
#    TC=TubeCollection()
#    TC.extend(d)
#    TC.update()
#    print TC.get_Nodes,d,TC,TC[0].TubeFcn
#    return TC

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
        cdef _Tube Tube
        for Tube in self:
            if Tube.key1 == key or Tube.key2 == key:
                return Tube
        raise KeyError
    
#    def __reduce__(self):
#        d=self[:]
#        print d
#        return rebuildTubeCollection,(d,)
    
cdef class _ControlVolume:
    pass

cdef class _Tube:
    pass

cpdef list collect_State_h(list CVList):
    cdef _ControlVolume CV
    return [CV.State.get_h() for CV in CVList]