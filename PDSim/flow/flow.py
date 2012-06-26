import numpy as np,copy
from PDSim.misc._listmath import listm
from _flow import _FlowPath,_FlowPathCollection
from _sumterms import sumterms_helper

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
    
    def calculate(self, Core):
        """
        Core is the main model core, it contains information that 
        is needed for the flow models
        """
        exists_keys=Core.CVs.exists_keys
        Tubes_Nodes=Core.Tubes.Nodes
        for FlowPath in self:
            ## Update the pointers to the states for the ends of the flow path
            if FlowPath.key1 in exists_keys:
                FlowPath.State1=Core.CVs[FlowPath.key1].State
            elif FlowPath.key1 in Tubes_Nodes:
                FlowPath.State1=Tubes_Nodes[FlowPath.key1]
            else:
                FlowPath.mdot=0.0
                #Doesn't exist, go to next flow
                continue                    
            
            if FlowPath.key2 in exists_keys:
                FlowPath.State2=Core.CVs[FlowPath.key2].State
            elif FlowPath.key2 in Tubes_Nodes:
                FlowPath.State2=Tubes_Nodes[FlowPath.key2]
            else:
                FlowPath.mdot=0.0
                #Doesn't exist, go to next flow
                continue
          
            #Calculate using the unpacked keyword arguments
            FlowPath.calculate(**FlowPath.MdotFcn_kwargs)
            
    def sumterms(self,Core):
        
        #Cached versions of the CV existence information
        exists_keys=Core.CVs.exists_keys #A list of the keys corresponding to the CV that exist
        omega=Core.omega
        
        #Call the Cython-version of the summer
        summerdm,summerdT = sumterms_helper(self,exists_keys,omega)
        
        return listm(summerdT),listm(summerdm)

def rebuildFlowPath(d):
    FP = FlowPath()
    for item in d:
        setattr(FP,item,d[item])
    return FP
    
class FlowPath(_FlowPath):
    """
    This is the class that contains the data and handlers for
    the flow models.  
    """
        
    def __init__(self,key1='',key2='',MdotFcn='',MdotFcn_kwargs={}):
        self.key1=key1
        self.key2=key2
        #Add the bound method in a wrapper that will pickle properly
        self.MdotFcn=MdotFcn
        self.MdotFcn_kwargs=MdotFcn_kwargs
    def update(self,d):
        if isinstance(d,dict):
            for k,v in d.iteritems():
                setattr(self,k,v)
    def calculate(self):
        """
        calculate
        """
        # Call the Cython routine in _Flow.pyx to evaluate the states
        # in preparation for the mass flow function
        self._calculate()
        
        # Pass off to the calculation function
        #
        # This function is probably in another class, so pass an 
        # instance of the current FlowPath class to the function
        #
        #Some other keyword arguments can be passed along to the function
        self.mdot=self.MdotFcn(self,**self.MdotFcn_kwargs)
    
    def __reduce__(self):
        return rebuildFlowPath,(self.__getstate__(),)
        
    def __getstate__(self):
        d={}
        d.update(self.__dict__.copy())
        d.update(self.__cdict__().copy())
        return d
        
    def __deepcopy__(self,memo):
        newFM=FlowPath()
        newdict = copy.copy(self.__dict__)
        newFM.__dict__.update(newdict)
        newFM.mdot=copy.copy(self.mdot)
        newFM.key_up=copy.copy(self.key_up)
        newFM.key_down=copy.copy(self.key_down)
        newFM.key1=copy.copy(self.key1)
        newFM.key2=copy.copy(self.key2)
        newFM.h_up=copy.copy(self.h_up)
        return newFM

if __name__=='__main__':
    FP=FlowPath()
    FL=FlowPathCollection()
    FL.append(FP)
    from copy import deepcopy
    f2=deepcopy(FL)
    print f2
    