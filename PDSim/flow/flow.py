import numpy as np,copy
from PDSim.misc._listmath import listm
from _flow import _FlowPathCollection, FlowPath
from _sumterms import sumterms_helper
from PDSim.flow.flow_models import FlowFunctionWrapper, PyFlowFunctionWrapper 
import cPickle

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
        return func.__get__(obj, cls)

import copy_reg
import types
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
            
class struct(object):
    def __init__(self,d):
        self.__dict__.update(d)

def rebuildFPC(d):
    """
    Used with cPickle to recreate the (empty) flow path collection class
    """
    FPC = FlowPathCollection()
    for Flow in d:
        FPC.append(Flow)
    return FPC
 
class FlowPathCollection(_FlowPathCollection):
    def __init__(self,d=None):
        if d is not None and isinstance(d,dict):
            self.__dict__.update(d)
    
    def __reduce__(self):
        return rebuildFPC,(self[:],)
    
    def get_deepcopy(self):
        """
        Using this method, the link to the mass flow function is broken
        """
        return [Flow.get_deepcopy() for Flow in self]
        FL=[]
        for Flow in self:
            FP=FlowPath()
            FP.update(Flow.__cdict__())
            FL.append(FP)
        return FL
    
    def deepcopy(self):
        newFPC = FlowPathCollection()
        newFPC.__dict__.update(copy.deepcopy(self.__dict__))
        return newFPC
        
    @property
    def N(self):
        return self.__len__()
    
        
    
            
    def sumterms(self,Core):

        #Call the Cython-version of the summer
        summerdm,summerdT = sumterms_helper(self,Core.CVs.exists_keys,Core.omega)
        
        return listm(summerdT),listm(summerdm)

    