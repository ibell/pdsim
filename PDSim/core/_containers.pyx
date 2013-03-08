from __future__ import division
cimport cython
    
cdef public enum STATE_VARS:
    STATE_VARS_TD
    STATE_VARS_TM
    
cdef class TubeCollection(list):
    
    def __init__(self):
        self._Nodes = {}
    
    cpdef dict get_Nodes(self):
        self.update()
        return self._Nodes
    
    cpdef arraym get_h(self):
        """
        Get an arraym instance with the enthalpies of each node in the Tubes
        collection.  In the same order as the indices of the enthalpies, but offset 
        by the number of control volumes in existence
        """
        return self.harray
        
    cpdef update_existence(self, int NCV):
        """
        Set the indices for each tube node in the array of enthalpies
        
        First index is equal to NCV since python (& c++) are 0-based indexing
        """
        cdef int i = NCV
        h = []
        for Tube in self:
            h.append(Tube.State1.h)
            h.append(Tube.State2.h)
            Tube.i1 = i
            Tube.i2 = i+1 
            i += 2
        self.harray = arraym(h)
        
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
    
    def __init__(self, int N):
        self.array_list = ['T','p','h','rho','V','dV','cp','cv','m','v',
                           'dpdT_constV','Q','xL','dudxL','drhodtheta', 
                           'dTdtheta', 'dmdtheta', 'dxLdtheta', 'summerdm', 
                           'summerdT', 'summerdxL', 'property_derivs']
        self.build_all(N)
        
    cdef build_all(self, int N):
        """
        Allocate the arrays, each the length of the number of CV in existence
        """
        for array_name in self.array_list:
            arr = arraym()
            arr.set_size(N)
            setattr(self,array_name,arr)
        
    cdef free_all(self):
        """
        Free all the arrays allocated
        """
        for array_name in self.array_list:
            if hasattr(self,array_name):
                delattr(self,array_name)
        
    cpdef update_size(self, int N):
        self.free_all()
        self.build_all(N)
         
    cpdef just_volumes(self, CVs, double theta):
        """
        Just calculate the volumes
        for each control volume.
        
        Parameters
        ----------
        CVs : list of control volumes
        theta : double
            Crank angle [radians]
        """
        
        cdef int N = len(CVs)
        
        for iCV in range(N):
            # Early-bind the control volume for speed
            CV = CVs[iCV]
            
            # Calculate the volume and derivative of volume - does not depend on
            # any of the other state variables
            self.V.data[iCV], self.dV.data[iCV] = CV.V_dV(theta, **CV.V_dV_kwargs)
    
    @cython.cdivision(True)    
    cpdef properties_and_volumes(self, CVs, double theta, int state_vars, arraym x):
        """
        Calculate all the required thermodynamic properties as well as the volumes
        for each control volume.
        
        Parameters
        ----------
        CVs : list of control volumes
        theta : double
            Crank angle [radians]
        state_vars : int
            Flag for the set of input variables - one of STATE_VARS_TM or STATE_VARS_TD defined in this module
        x : arraym
            List of state variables corresponding to the state_vars flag
        """
        cdef StateClass State
        cdef int N = len(CVs)
        cdef int iCV, iVar
        cdef arraym T,rho,m
        cdef arraym var1, var2 #Arrays to hold the state values for T,rho for instance
        
        #Calculate the volumes
        self.just_volumes(CVs,theta)
        
        # Split the state variable array into chunks
        self.T = x.slice(0,N)
        if state_vars == STATE_VARS_TM:
            m = x.slice(N, 2*N)
            self.m.set_data(m.data,N)
            for iCV in range(N):
                self.rho.data[iCV] = self.m.data[iCV]/self.V.data[iCV]
        elif state_vars == STATE_VARS_TD:
            rho = x.slice(N, 2*N)
            self.rho.set_data(rho.data, N)
            for iCV in range(N):
                self.m.data[iCV] = self.rho.data[iCV]*self.V.data[iCV]
        
        for iCV in range(N):
            self.v.data[iCV] = 1/self.rho.data[iCV]
        
        for iCV in range(N):
            # Early-bind the control volume and State for speed
            CV = CVs[iCV]
            State = CV.State

            # Update the CV state variables using temperature and density
            State.update_Trho(self.T.data[iCV], self.rho.data[iCV])
            
            self.p.data[iCV] = State.get_p()
            self.h.data[iCV] = State.get_h()
            self.cp.data[iCV] = State.get_cp()
            self.cv.data[iCV] = State.get_cv()
            self.dpdT_constV.data[iCV] = State.get_dpdT()
        
        self.N = N
        self.state_vars = state_vars
    
    cpdef calculate_flows(self, FlowPathCollection Flows, arraym harray, Core):
        """
        Calculate the flows between tubes and control volumes and sum up the 
        flow-related terms
        
        Parameters
        ----------
        Flows: :class:`PDSim.core.flow._flow.FlowPathCollection` instance
        harray: :class:`PDSim.misc.datatypes.arraym` instance
            An array,
        
        """
        cdef int i
        cdef arraym summerdm, summerdT
        
        Flows.calculate(harray)
        self.summerdT, self.summerdm = Flows.sumterms(Core)
    
    @cython.cdivision(True)
    cpdef calculate_derivs(self, double omega, bint has_liquid):
        
        cdef double m,T,cv,xL,dV,V,v,summerdxL,summerdm,summerdT
        
        self.omega = omega
        
        #Set some variables for the oil-flooded case which is not yet supported
        self.xL = arraym()
        self.xL.set_size(self.N)
        self.dudxL = arraym()
        self.dudxL.set_size(self.N)
        self.summerdxL = arraym()
        self.summerdxL.set_size(self.N)
        
        #The derivative arrays
        self.dxLdtheta = arraym()
        self.dxLdtheta.set_size(self.N)
        self.dTdtheta = arraym()
        self.dTdtheta.set_size(self.N)
        self.drhodtheta = arraym()
        self.drhodtheta.set_size(self.N)
        
        #Actually calculate the derivatives
        self.dmdtheta = self.summerdm
        
        #Loop over the control volumes
        for i in range(self.N):
            #For compactness, pull the data from the arrays
            m = self.m.data[i]
            h = self.h.data[i]
            T = self.T.data[i]
            Q = self.Q.data[i]
            omega = self.omega
            rho = self.rho.data[i]
            cv = self.cv.data[i]
            xL = self.xL.data[i]
            dV = self.dV.data[i]
            V = self.V.data[i]
            v = self.v.data[i]
            dudxL = self.dudxL.data[i]
            dpdT = self.dpdT_constV.data[i]
            summerdxL = self.summerdxL.data[i]
            summerdm = self.summerdm.data[i]
            summerdT = self.summerdT.data[i]
            dmdtheta = self.dmdtheta.data[i]
            
            self.dxLdtheta.data[i] = 1.0/m*(summerdxL-xL*dmdtheta);    dxLdtheta = self.dxLdtheta.data[i]            
            self.dTdtheta.data[i] = 1.0/(m*cv)*(-1.0*T*dpdT*(dV-v*dmdtheta)-m*dudxL*dxLdtheta-h*dmdtheta+Q/omega+summerdT)
            self.drhodtheta.data[i] = 1.0/V*(dmdtheta-rho*dV)
        
        # Create the array of output values
        self.property_derivs = self.dTdtheta.copy()
        if self.state_vars == STATE_VARS_TM:
            self.property_derivs.extend(self.dmdtheta)
        elif self.state_vars == STATE_VARS_TD:
            self.property_derivs.extend(self.drhodtheta)
            
    cpdef copy(self):
        CVA = CVArrays(self.T.N)
        #Loop over the names of the arrays
        for array_name in self.array_list:
            #Get the array from this class
            arr = getattr(self,array_name)
            #Put a copy of it into the new class
            setattr(CVA,array_name,<arraym>arr.copy())
        
        return CVA
        
cdef class ControlVolume(object):
    """
    This is a class that contains all the code for a given control volume.  
    
    It includes the code for calculation of volumes and others.
    
    Parameters
    ----------
    key : str
        The string of the key for this control volume
    VdVFcn : function, (future: VolumeFunction class)
    initialState : :class:`State <CoolProp.State.State>` instance
    exists : bool
        ``True`` if control volume exists, ``False`` otherwise
    VdVFcn_kwargs : dict
        Keyword arguments that can be passed to the VdVFcn 
    discharge_becomes : str
        The key of the chamber that this control volume becomes at the 
        discharge angle (scroll compressor only)
    becomes : str or list
        The key of the control volume that this CV becomes in the next revolution,
        or a list of keys of control volumes that take on the values of this
        CV
    """
    
    def __init__(self, 
                 str key, 
                 object VdVFcn, 
                 StateClass initialState, 
                 bint exists=True,
                 dict VdVFcn_kwargs={}, 
                 str discharge_becomes=None, 
                 object becomes = None):

        self.key = key.encode('ascii')
        self.V_dV = VdVFcn
        self.State = initialState
        self.exists = exists
        self.V_dV_kwargs = VdVFcn_kwargs #Keyword-arguments that can get passed to volume function
        self.discharge_becomes = discharge_becomes.encode('ascii') if discharge_becomes is not None else key.encode('ascii')
        self.becomes = becomes if becomes is not None else key.encode('ascii')
    
    def __reduce__(self):
        #TODO: fix me
        return rebuildCV,(self.__getstate__().copy(),)
    
    def __getstate__(self):
        #TODO: fix me
        d=self.__dict__
        d['State']=self.State
        return d.copy()
    
    def __setstate__(self, d):
        #TODO: fix me
        for item in d:
            setattr(self,item,d[item])
        
    def __deepcopy__(self):
        #TODO: fix me
        import copy
        return copy.deepcopy(self)
    
def rebuildCV(dict d):
    CV = ControlVolume(d.pop('key'),d.pop('V_dV'),d.pop('State'))
    for item in d:
        setattr(CV,item,d[item])
    return CV
        
cdef class ControlVolumeCollection(object):
    """
    ControlVolumeCollection is class to hold all the control volumes
    """
    def __init__(self):
        self.keys = []
        self.CVs = []
        
    def __reduce__(self):
        #TODO: rewrite me
        return rebuildCVCollection,(self.__getstate__(),)
    
    def __getstate__(self):
        #TODO: rewrite me
        import copy
        CVs = [copy.copy(item) for k,item in self.iteritems()]
        return CVs

    def __setstate__(self, CVs):
        #TODO: rewrite me
        for CV in CVs:
            self[CV.key]=CV
            
    def __getitem__(self, k):
        """
        Can index based on integer index or string key
        """
        if k in self.keys:
            return self.CVs[self.keys.index(k)]
        elif k in range(self.N):
            return self.CVs[k]
        else:
            raise KeyError('Your key [{key:s}] is invalid'.format(key = k))
            
    cpdef add(self, ControlVolume CV):
        """
        Add a control volume to the list of control volumes
        
        Parameters
        ----------
        CV : :class:`ControlVolume <PDSim.core.containers.ControlVolume>' instance
        
        """
        if CV.key in self.keys:
            raise ValueError('Your CV key [{key:s}] is already in use'.format(CV.key))
        else:
            self.CVs.append(CV)
            self.keys.append(CV.key)
    
    cpdef rebuild_exists(self):
        
        # For all CV - whether they exist or not
        # both indices and keys are in the same order
        self.indices = range(0,len(self.keys))
        
        self.exists_indices = [i for i in self.indices if self.CVs[i].exists]
        self.exists_keys = [self.keys[i] for i in self.exists_indices]
        self.exists_CV = [self.CVs[i] for i in self.exists_indices]
        if len(self.exists_keys) == 0:
            return
        
        self.Nodes = dict([(CV.key, CV.State) for CV in self.exists_CV])
        self.N = len(self.CVs)
        self.Nexist = len(self.exists_CV)
    
    def index(self,key):
        return self.keys.index(key)
    
    @property
    def T(self):
        """
        Temperature for each CV that exists
        """
        cdef ControlVolume CV
        return [CV.State.get_T() for CV in self.exists_CV]
    
    @property
    def p(self):
        """
        Pressure for each CV that exists
        """
        cdef ControlVolume CV
        return [CV.State.get_p() for CV in self.exists_CV]
    
    @property
    def rho(self):
        """
        Density for each CV that exists
        """
        cdef ControlVolume CV
        return [CV.State.get_rho() for CV in self.exists_CV]
    
    @property
    def h(self):
        """
        Enthalpy for each CV that exists
        """
        cdef ControlVolume CV
        return [CV.State.get_h() for CV in self.exists_CV]
    
    @property
    def cp(self):
        """
        Specific heat at constant volume for each CV that exists
        """
        cdef ControlVolume CV
        return [CV.State.get_cp() for CV in self.exists_CV]
    
    @property
    def cv(self):
        """
        Specific heat at constant volume for each CV that exists
        """
        cdef ControlVolume CV
        return [CV.State.get_cv() for CV in self.exists_CV]
    
    @property
    def dpdT(self):
        """
        Derivative of pressure with respect to temperature at constant volume for each CV that exists
        """
        cdef ControlVolume CV
        return [CV.State.get_dpdT() for CV in self.exists_CV]
    
    cpdef updateStates(self, str name1, arraym array1, str name2, arraym array2):
#        if not len(array1) == len(array2) or not len(array2)==len(self.exists_CV):
#            raise AttributeError('length of arrays must be the same and equal number of CV in existence')
        keys = self.exists_keys
        # Update each of the states of the control volume
        for CV,v1,v2 in zip(self.exists_CV, array1, array2):
            CV.State.update({name1:v1,name2:v2})
     
    cpdef volumes(self, double theta, bint as_dict = False):
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
            
        #Loop over the control volumes that exist 
        V_dV=[CV.V_dV(theta, **CV.V_dV_kwargs) for CV in self.exists_CV] #See below for Vfunc
        V,dV=zip(*V_dV)
        if not as_dict:
            return arraym(V),arraym(dV)
        else:
            V_dict = {key:_V for key,_V in zip(self.exists_keys,V)}
            dV_dict = {key:_dV for key,_dV in zip(self.exists_keys,dV)}
            return V_dict, dV_dict

def rebuildCVCollection(CVs):
    CVC = ControlVolumeCollection()
    for CV in CVs:
        CVC[CV.key]=CV
    return CVC

cdef class _Tube:
    pass

cpdef list collect_State_h(list CVList):
    cdef _ControlVolume CV
    return [CV.State.get_h() for CV in CVList]