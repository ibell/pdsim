
from __future__ import division
cimport cython

cdef public enum STATE_VARS:
    STATE_VARS_TD
    STATE_VARS_TM
    
cdef class Tube(object):
    """
    A tube is a component of the model that allows for heat transfer and pressure drop.
    
    With this class, the state of at least one of the points is fixed.  For instance, at the inlet of the compressor, the state well upstream is quasi-steady.
    """
    def __init__(self,key1,key2,L,ID,State1=None,State2=None,OD=-1,fixed=-1,TubeFcn=None,mdot=-1,exists=True):
        """
        
        Parameters
        ----------
        key1 : str
            Key for the upstream flow node
        key2 : str
            Key for the downstream flow node
        L : float
            Length of the tube [m] 
        ID : float
            Internal diameter of the tube [m]
        State1 : :class:`State <CoolProp.State.State>` instance
            Upstream state
        State2 : :class:`State <CoolProp.State.State>` instance
            Downstream state
        OD : float
            Outer diameter of the tube [m]
        fixed : int
            Which one of the node is fixed, one of ``1`` or ``2``
        TubeFun : function
            A function that will be called for the tube
        mdot : float
            The mass flow rate [kg/s]
        exists : boolean
            ``True`` if the tube exists
        """
        self.key1 = key1
        self.key2 = key2
        self.fixed = fixed
        
        #: Additional heat to be added to the tube
        self.Q_add = 0.0
        
        #: Fixed heat transfer coefficient if desired (if less than zero will use correlation - default)
        self.alpha = -1.0
        
        self.exists = exists
        if fixed<0:
            raise AttributeError("You must provide an integer value for fixed, either 1 for Node 1 fixed, or 2 for Node 2 fixed.")
        if fixed==1 and isinstance(State1,StateClass) and State2==None:
            #Everything good
            self.State1=State1
            self.State2=State1.copy()
        elif fixed==2 and isinstance(State2,StateClass) and State1==None:
            #Everything good
            self.State2=State2
            self.State1=State2.copy()
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
    
    cpdef arraym get_p(self):
        """
        Get an arraym instance with the pressures of each node in the Tubes
        collection.  In the same order as the indices of the pressures, but offset 
        by the number of control volumes in existence
        """
        return self.parray
    
    cpdef arraym get_T(self):
        """
        Get an arraym instance with the enthalpies of each node in the Tubes
        collection.  In the same order as the indices of the enthalpies, but offset 
        by the number of control volumes in existence
        """
        return self.Tarray
        
    cpdef update_existence(self, int NCV):
        """
        Set the indices for each tube node in the array of enthalpies
        
        First index is equal to NCV since python (& c++) are 0-based indexing
        """
        cdef int i = NCV
        h,p,T = [],[],[]
        for Tube in self:
            h.append(Tube.State1.h)
            h.append(Tube.State2.h)
            p.append(Tube.State1.p)
            p.append(Tube.State2.p)
            T.append(Tube.State1.T)
            T.append(Tube.State2.T)
            Tube.i1 = i
            Tube.i2 = i+1 
            i += 2
        self.harray = arraym(h)
        self.parray = arraym(p)
        self.Tarray = arraym(T)
        
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

cdef class CVScore(object):
    """
    The base class for all control volumes

    In the derived class, before anything is done, you must set the parameter array_list as a list of strings, each
    entry in the list should be the name of an arraym instance that will be stored in the class
    """

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
        """
        Update the size of the arraym instances in this class
        """
        self.free_all()
        self.build_all(N)

    cpdef copy(self):
        """
        Make copies of all of the arraym instances in this class
        """
        CVA = type(self)(self.T.N) #type(self) gets the derived class' type
        # Loop over the names of the arrays
        for array_name in self.array_list:
            # Get the array from this class
            arr = getattr(self,array_name)
            # Put a copy of it into the new class
            setattr(CVA, array_name, <arraym>arr.copy())
        return CVA

    cpdef calculate_flows(self, FlowPathCollection Flows, arraym harray, arraym parray, arraym Tarray):
        """
        Calculate the flows between tubes and control volumes and sum up the
        flow-related terms

        Loads the arraym instances ``summerdT`` and ``summerdm`` of this class

        These terms are defined by

        .. math::

            \\mathrm{summerdm} = \\sum  \\frac{\\dot m}{\\omega}

        and

        .. math::

            \\mathrm{summerdT} = \\sum  \\frac{\\dot m h}{\\omega}

        where the signs are dependent on whether the flow is into or out of the
        given control volume

        Parameters
        ----------
        Flows : :class:`FlowPathCollection <PDSim.flow.flow.FlowPathCollection>` instance
        harray : :class:`arraym <PDSim.misc.datatypes.arraym>` instance
        parray : :class:`arraym <PDSim.misc.datatypes.arraym>` instance
        Tarray : :class:`arraym <PDSim.misc.datatypes.arraym>` instance
        """

        Flows.calculate(harray, parray, Tarray)
        Flows.sumterms(self.summerdT, self.summerdm)

cdef class CVArrays(CVScore):
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
         
    cpdef just_volumes(self, list CVs, double theta):
        """
        Just calculate the volumes
        for each control volume.
        
        Parameters
        ----------
        CVs : list of control volumes
        theta : double
            Crank angle [radians]
        """

        cdef int N = len(CVs), iCV
        
        for iCV in range(N):
            # Early-bind the control volume for speed
            CV = <ControlVolume>(CVs[iCV])
            
            # Calculate the volume and derivative of volume - does not depend on
            # any of the other state variables
            self.V.data[iCV], self.dV.data[iCV] = CV.V_dV(theta, **CV.V_dV_kwargs)
    
    @cython.cdivision(True)    
    cpdef properties_and_volumes(self, list CVs, double theta, int state_vars, arraym x):
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
        x : :class:`arraym <PDSim.misc.datatypes.arraym>` instance
            List of state variables corresponding to the state_vars flag
        """
        cdef StateClass State
        cdef int N = len(CVs)
        cdef int iCV, iVar, i
        
        #Calculate the volumes
        self.just_volumes(CVs,theta)
        
        # Split the state variable array into chunks
        for i in range(N):
            self.T.data[i] = x.data[i]
            
        if state_vars == STATE_VARS_TM:
            i = 0
            for j in xrange(N, 2*N):
                self.m.data[i] = x.data[j]
                i += 1
            for iCV in range(N):
                self.rho.data[iCV] = self.m.data[iCV]/self.V.data[iCV]
        elif state_vars == STATE_VARS_TD:
            i = 0
            for j in xrange(N, 2*N):
                self.rho.data[i] = x.data[j]
                i += 1
            for iCV in range(N):
                self.m.data[iCV] = self.rho.data[iCV]*self.V.data[iCV]
        
        for iCV in range(N):
            self.v.data[iCV] = 1/self.rho.data[iCV]
        
        for iCV in range(N):
            # Early-bind the State for speed
            State = (<ControlVolume>(CVs[iCV])).State

            # Update the CV state variables using temperature and density
            State.update_Trho(self.T.data[iCV], self.rho.data[iCV])
            
            self.p.data[iCV] = State.get_p()
            self.h.data[iCV] = State.get_h()
            self.cp.data[iCV] = State.get_cp()
            self.cv.data[iCV] = State.get_cv()
            self.dpdT_constV.data[iCV] = State.get_dpdT()
        
        self.N = N
        self.state_vars = state_vars
    
    @cython.cdivision(True)
    cpdef calculate_derivs(self, double omega, bint has_liquid):
        
        cdef double m,T,cv,xL,dV,V,v,summerdxL,summerdm,summerdT
        cdef int i
        
        self.omega = omega
        
        #Actually calculate the derivatives
        self.dmdtheta = self.summerdm
        self.property_derivs = arraym()
        self.property_derivs.set_size(self.N*2)
        
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
            
            self.dxLdtheta.data[i] = 1.0/m*(summerdxL-xL*dmdtheta)    
            dxLdtheta = self.dxLdtheta.data[i]            
            self.dTdtheta.data[i] = 1.0/(m*cv)*(-1.0*T*dpdT*(dV-v*dmdtheta)-m*dudxL*dxLdtheta-h*dmdtheta+Q/omega+summerdT)
            self.drhodtheta.data[i] = 1.0/V*(dmdtheta-rho*dV)
        
        #  Create the array of output values
        for i in range(self.N):
            self.property_derivs.set_index(i, self.dTdtheta.data[i])
            
            if self.state_vars == STATE_VARS_TM:
                self.property_derivs.set_index(i + self.N, self.dmdtheta.data[i])
            elif self.state_vars == STATE_VARS_TD:
                self.property_derivs.set_index(i + self.N, self.drhodtheta.data[i])
        
cdef class ControlVolume(object):
    """
    This is a class that contains all the code for a given control volume.  
    
    It includes the code for calculation of volumes and others.
    """
    
    def __init__(self, 
                 str key, 
                 object VdVFcn, 
                 StateClass initialState, 
                 bint exists = True,
                 dict VdVFcn_kwargs = {},
                 str discharge_becomes = None,
                 object becomes = None):
        """
        Parameters
        ----------
        key : str
            The string of the key for this control volume
        VdVFcn : function, (future: VolumeFunction class)
        initialState : :class:`State <CoolProp.State.State>` instance
        exists : bool
            ``True`` if control volume exists, ``False`` otherwise
        VdVFcn_kwargs : dict
            Dictionary of keyword arguments that can be passed to the VdVFcn
        discharge_becomes : str
            The key of the chamber that this control volume becomes at the 
            discharge angle (scroll compressor only)
        becomes : str or list
            The key of the control volume that this CV becomes in the next revolution,
            or a list of keys of control volumes that take on the values of this
            CV
        """

        self.key = key.encode('ascii')
        self.V_dV = VdVFcn
        self.State = initialState
        self.exists = exists
        self.V_dV_kwargs = VdVFcn_kwargs
        self.discharge_becomes = discharge_becomes.encode('ascii') if discharge_becomes is not None else key.encode('ascii')
        self.becomes = becomes if becomes is not None else key.encode('ascii')
        
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
        import copy
        CVs = [copy.copy(CV) for CV in self.CVs]
        return CVs

    def __setstate__(self, CVs):
        #TODO: rewrite me
        for CV in CVs:
            self[CV.key]=CV
            
    cpdef at(self, int i):
        """
        Return the control volume at the given index
        """
        return self.CVs[i]
    
    def __getitem__(self, k):
        """
        Can index based on integer index or string key
        """
        if k in self.keys:
            return self.CVs[self.keys.index(k)]
        elif k in range(self.N):
            return self.CVs[k]
        else:
            raise KeyError('Your key [{key:s}] of type [{_type:s}] is invalid'.format(key = k,_type = str(type(k))))
            
    def __len__(self):
        return len(self.CVs)
            
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
        """
        Rebuild all the internal lists that hold the indices, keys, and control volumes
        """
        
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
        CVC.add(CV)
    CVC.rebuild_exists()
    return CVC

cpdef list collect_State_h(list CVList):
    cdef ControlVolume CV
    return [CV.State.get_h() for CV in CVList]