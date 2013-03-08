import copy
import collections
from cPickle import dumps,loads
from CoolProp.State import State

from _containers import TubeCollection, CVArrays
from PDSim.misc.datatypes import arraym
    
class Tube(object):
    """
    A tube is a component of the model that allows for heat transfer and pressure drop.
    
    With this class, the state of at least one of the points is fixed.  For instance, at the inlet of the compressor, the state well upstream is quasi-steady.
    """
    def __init__(self,key1,key2,L,ID,State1=None,State2=None,OD=-1,fixed=-1,TubeFcn=None,mdot=-1,exists=True):
        self.key1 = key1
        self.key2 = key2
        self.fixed = fixed
        
        #: Additional heat to be added to the tube
        self.Q_add = 0.0
        
        #: Fixed heat transfer coefficient if desired (if less than zero will use correlation - default)
        self.alpha = -1.0
        
        self.exists = exists
        if fixed<0:
            raise AttributeError(textwrap.dedent("""You must provide an integer 
            value for fixed, either 1 for Node 1 fixed, or 2 for Node 2 fixed.  
            You provided None (or didn\'t include the parameter"""))
        if fixed==1 and isinstance(State1,State) and State2==None:
            #Everything good
            self.State1=State1
            self.State2=State(self.State1.Fluid,{'T':self.State1.T,'D':self.State1.rho})
        elif fixed==2 and isinstance(State2,State) and State1==None:
            #Everything good
            self.State2=State2
            self.State1=State(self.State2.Fluid,{'T':self.State2.T,'D':self.State2.rho})
        else:
            raise AttributeError('Incompatibility between the value for fixed and the states provided')
            
        self.TubeFcn=TubeFcn
        if mdot<0:
            self.mdot=0.010
            print('Warning: mdot not provided to Tube class constructor, guess value of '+str(self.mdot)+' kg/s used')
        else:
            self.mdot=mdot
        self.L=L
        self.ID=ID
        self.OD=OD
        
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
        """
        A list of all the nodes associated with the control volumes
        """
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
        """
        Temperature for each CV that exists
        """
        return [CV.State.get_T() for CV in self._exists_CV]
    
    @property
    def p(self):
        """
        Pressure for each CV that exists
        """
        return [CV.State.get_p() for CV in self._exists_CV]
    
    @property
    def rho(self):
        """
        Density for each CV that exists
        """
        return [CV.State.get_rho() for CV in self._exists_CV]
    
    @property
    def h(self):
        """
        Enthalpy for each CV that exists
        """
        return [CV.State.get_h() for CV in self._exists_CV]
    
    @property
    def cp(self):
        """
        Specific heat at constant volume for each CV that exists
        """
        return [CV.State.get_cp() for CV in self._exists_CV]
    
    @property
    def cv(self):
        """
        Specific heat at constant volume for each CV that exists
        """
        return [CV.State.get_cv() for CV in self._exists_CV]
    
    @property
    def dpdT(self):
        """
        Derivative of pressure with respect to temperature at constant volume for each CV that exists
        """
        return [CV.State.get_dpdT() for CV in self._exists_CV]
    ## ---- End property callbacks --------
    
    def updateStates(self,name1,array1,name2,array2,keys=None):
#        if not len(array1) == len(array2) or not len(array2)==len(self.exists_CV):
#            raise AttributeError('length of arrays must be the same and equal number of CV in existence')
        if keys is None:
            keys=self.exists_keys
        # Update each of the states of the control volume
        for CV,v1,v2 in zip(self._exists_CV, array1, array2):
            CV.State.update({name1:v1,name2:v2})
     
    def volumes(self,theta, as_dict = False):
        """
        Each control volume class must define a function V_dV (through a pointer) 
        that defines the volume and derivative of volume with respect to the 
        independent variable.  The function that V_dV points to MUST be of the form
        
        V,dV=V_dV(theta,**kwargs)
        
        If the parameter V_dV_kwargs is passed to the class constructor, these keyword 
        arguments will be unpacked into the volume function call.  Useful for passing 
        a flag to a given function
        
        Parameters
        ----------
        as_dict : boolean, optional
            If ``True``, return the volumes and derivatives of volumes as a dictionary
            
        Returns
        -------
        A tuple of volumes and derivatives of volumes as arraym instances
        
        """
            
        def func(CV):
            return CV.V_dV(theta,**CV.V_dV_kwargs)
        #Loop over the control volumes that exist 
        V_dV=map(func,self.exists_CV)
        V,dV=zip(*V_dV)
        if not as_dict:
            return arraym(V),arraym(dV)
        else:
            V_dict = {key:_V for key,_V in zip(self.exists_keys,V)}
            dV_dict = {key:_dV for key,_dV in zip(self.exists_keys,dV)}
            return V_dict, dV_dict

def rebuildCV(d):
    CV = ControlVolume(d.pop('key'),d.pop('V_dV'),d.pop('State'))
    for item in d:
        setattr(CV,item,d[item])
    return CV

class ControlVolume(object):
    """
    This is a class that contains all the code for a given control volume.  
    
    It includes the code for calculation of volumes and others.
    """
    
    def __init__(self, key, VdVFcn, initialState, exists=True,
                 VdVFcn_kwargs={}, discharge_becomes=None, becomes=None):
        #_ControlVolume.__init__(self)
        self.State=initialState
        self.exists=exists
        self.key=key
        
        self.V_dV = VdVFcn
        self.V_dV_kwargs = VdVFcn_kwargs #Keyword-arguments that can get passed to volume function
        self.discharge_becomes = discharge_becomes if discharge_becomes is not None else key
        self.becomes=becomes if becomes is not None else key
    
    def __reduce__(self):
        return rebuildCV,(self.__getstate__().copy(),)  
    
    def __getstate__(self):
        d=self.__dict__
        d['State']=self.State
        return d.copy()
    
    def __setstate__(self, d):
        for item in d:
            setattr(self,item,d[item])
        
    def __deepcopy__(self):
        return copy.deepcopy(self)