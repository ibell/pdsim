
cimport cython

from CoolProp.State import State as StateClass
from CoolProp.State cimport State as StateClass
    
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
            Flag for the set of input variables - one of STATE_VARS_TM or STATE_VARS_TD
        x : arraym
            List of state variables corresponding to the state_vars flag
        """
        cdef StateClass State
        cdef tuple V_dV
        cdef int N = len(CVs)
        cdef int iCV, iVar
        cdef arraym T,rho,m
        cdef arraym var1, var2 #Arrays to hold the state values for T,rho for instance
        
        # Split the state variable array into chunks
        T = x.slice(0,N)
        if state_vars == STATE_VARS_TM:
            m = x.slice(N, 2*N)
        elif state_vars == STATE_VARS_TD:
            rho = x.slice(N, 2*N)
        
        for iCV in range(N):
            # Early-bind the control volume and State for speed
            CV = CVs[iCV]
            State = CV.State
            
            # Calculate the volume and derivative of volume - does not depend on
            # any of the other state variables
            self.V.data[iCV], self.dV.data[iCV] = CV.V_dV(theta, **CV.V_dV_kwargs)
            
            # Update the state variables
            if state_vars == STATE_VARS_TM:
                # Update the CV state variables using temperature and mass
                State.update_Trho(T.data[iCV], m.data[iCV] / self.V.data[iCV])
            elif state_vars == STATE_VARS_TD:
                # Update the CV state variables using temperature and density
                State.update_Trho(T.data[iCV], rho.data[iCV])
            
            self.T.data[iCV] = State.get_T()
            self.p.data[iCV] = State.get_p()
            self.rho.data[iCV] = State.get_rho()
            self.h.data[iCV] = State.get_h()
            self.cp.data[iCV] = State.get_cp()
            self.cv.data[iCV] = State.get_cv()
            self.dpdT_constV.data[iCV] = State.get_dpdT()
            
            self.m.data[iCV] = self.rho.data[iCV] * self.V.data[iCV]
            self.v.data[iCV] = 1 / self.rho.data[iCV]
        
        self.N = N
        self.state_vars = state_vars
    
    cpdef calculate_flows(self, Flows, harray, Core):
        """
        Calculate the flows between tubes and control volumes and sum up the flow-related terms 
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
        
    
cdef class _ControlVolume:
    pass

cdef class _Tube:
    pass

cpdef list collect_State_h(list CVList):
    cdef _ControlVolume CV
    return [CV.State.get_h() for CV in CVList]