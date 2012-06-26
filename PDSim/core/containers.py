from PDSim.misc._listmath import listm
import copy
import collections
from cPickle import dumps,loads

def rebuildCVCollection(CVs):
    CVC = ControlVolumeCollection()
    for CV in CVs:
        CVC[CV.key]=CV
    return CVC
    
class ControlVolumeCollection(collections.OrderedDict):
    """
    ControlVolumeCollection is an extended dictionary with some PDSim related functions added
    """
    def __init__(self):
        collections.OrderedDict.__init__(self)
        #delattr(collections.OrderedDict,'__reduce__')
        
    def __reduce__(self):
        return rebuildCVCollection,(self.__getstate__(),)
    
    def __getstate__(self):
        CVs = [copy.copy(item) for k,item in self.iteritems()]
        return CVs

    def __setstate__(self, CVs):
        for CV in CVs:
            self[CV.key]=CV
    
    def rebuild_exists(self):
        # For all CV - whether they exist or not
        # both _indices and _keys are in the same order, thanks to the
        # use of the OrderedDict
        self._keys = self.keys()
        self._indices = range(0,len(self._keys))
        
        #For the CV in existence
        self._exists_keys = [k for k in self._keys if self[k].exists==True]
        if len(self._exists_keys)==0:
            return
        #Get the key for each control volume - sorted in the same order as _exists_indices
        self._exists_indices = [self._keys.index(k) for k in self._exists_keys]
        self._exists_CV = [self[k] for k in self._exists_keys]
        self._Nodes = dict([(CV.key, CV.State) for CV in self.exists_CV])
        self._Nexist = self.exists_CV.__len__()
    
    @property
    def Nodes(self):
        return self._Nodes
    
    def index(self,key):
        return self._keys.index(key)
    
    @property
    def exists_keys(self):
        return self._exists_keys
    
    @property
    def exists_indices(self):
        return self._exists_indices
    
    @property
    def N(self):
        return self.__len__()
    
    @property
    def Nexist(self):
        return self._Nexist
    
    @property
    def exists_CV(self):
        return self._exists_CV
            
    @property
    def T(self):
        return [CV.State.get_T() for CV in self._exists_CV]
    
    @property
    def p(self):
        return [CV.State.get_p() for CV in self._exists_CV]
    
    @property
    def rho(self):
        return [CV.State.get_rho() for CV in self._exists_CV]
    
    @property
    def h(self):
        #return collect_State_h(self._exists_CV)
#        hlist2= [CV.State.get_h() for CV in self._exists_CV]
#        print hlist1
#        print hlist2
        return [CV.State.get_h() for CV in self._exists_CV]
    
    @property
    def cp(self):
        return [CV.State.get_cp() for CV in self._exists_CV]
    
    @property
    def cv(self):
        return [CV.State.get_cv() for CV in self._exists_CV]
    
    @property
    def dpdT(self):
        return [self[k].State.dpdT for k in self._exists_keys]
    ## ---- End property callbacks --------
    
    def updateStates(self,name1,array1,name2,array2,keys=None):
#        if not len(array1) == len(array2) or not len(array2)==len(self.exists_CV):
#            raise AttributeError('length of arrays must be the same and equal number of CV in existence')
        if keys==None:
            keys=self.exists_keys
        # Update each of the states of the control volume
        for k,v1,v2 in zip(keys,array1,array2):
            self[k].State.update({name1:v1,name2:v2})
     
    def volumes(self,theta):
        """
        Each control volume class must define a function V_dV (through a pointer) 
        that defines the volume and derivative of volume with respect to the 
        independent variable.  The function that V_dV points to MUST be of the form
        
        V,dV=V_dV(theta,**kwargs)
        
        If the parameter V_dV_kwargs is passed to the class constructor, these keyword 
        arguments will be unpacked into the volume function call.  Useful for passing 
        a flag to a given function
        """
            
        def func(CV):
            return CV.V_dV(theta,**CV.V_dV_kwargs)
        #Loop over the control volumes that exist 
        V_dV=map(func,self.exists_CV)
        V,dV=zip(*V_dV)
        return listm(V),listm(dV)

def rebuildCV(d):
    CV = ControlVolume(d.pop('key'),d.pop('V_dV'),d.pop('State'))
    for item in d:
        setattr(CV,item,d[item])
    return CV
    
class ControlVolume(object):
    """
    This is a class that contains all the code for a given control volume.  
    
    It includes the code for calculation of volumes and others.  Some functions 
    must be overloaded on instantiation
    """
    
    def __init__(self,key,VdVFcn,initialState,exists=True,VdVFcn_kwargs={},discharge_becomes='',becomes=''):
        #_ControlVolume.__init__(self)
        self.State=initialState
        self.exists=exists
        self.key=key
        self.V_dV=VdVFcn
        self.V_dV_kwargs=VdVFcn_kwargs #Keyword-arguments that can get passed to volume function
        self.discharge_becomes=discharge_becomes
        self.becomes=becomes
    
    def __reduce__(self):
        return rebuildCV,(self.__getstate__().copy(),)
    
    def __getstate__(self):
        d=self.__dict__
        d['State']=self.State
        return d.copy()
#    
#    
    def __setstate__(self, d):
        for item in d:
            setattr(self,item,d[item])
#            
#    def __getinitargs__(self):
#        return (self.key, self.V_dV, self.State)
        
#    def __reduce__(self):
#        print 'reduce'
#        d=self.__dict__
#        d['State']=self.State
#        return rebuildCV, (d, )
        
    def __deepcopy__(self):
        return copy.deepcopy(self)