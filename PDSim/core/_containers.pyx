
cdef class _ControlVolume:
    pass
#    def __reduce__(self):
#        d={}
#        return rebuildCV,(d,)

#def rebuildCV(d):
#    """
#    Used with cPickle to recreate the (empty) control volume class
#    """
#    CV=_ControlVolume()
#    return CV

cpdef list collect_State_h(list CVList):
    cdef _ControlVolume CV
    return [CV.State.get_h() for CV in CVList]