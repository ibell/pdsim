
#Other imports in _flow.pxd
import copy
import cPickle
from flow_models import PyFlowFunctionWrapper

from CoolProp.State import State as StateClass
from CoolProp.State cimport State as StateClass

import cython

cpdef tuple sumterms_given_CV(bytes key, list Flows):
    """
    A function to sum all the mdot terms for a given control volume
    
    Searches the list of flows and for each element in the flows, checks whether
    the key matches one of upstream or downstream key
    
    Parameters
    ----------
    key: string
    Flows: FlowPathCollection instance 
    """
    cdef int i
    cdef FlowPath Flow
    cdef double summer_mdot = 0.0, summer_mdoth = 0.0
    
    for i in range(len(Flows)):
        Flow = <FlowPath>Flows[i]

        if not Flow.exists or abs(Flow.mdot)<1e-12:
            continue
                
        if Flow.key_down == key:
            summer_mdot+=Flow.mdot
            summer_mdoth+=Flow.mdot*Flow.h_up
        
        elif Flow.key_up == key:
            summer_mdot-=Flow.mdot
            summer_mdoth-=Flow.mdot*Flow.h_up
            
        else:
            continue
        
    return summer_mdot, summer_mdoth

class struct(object):
    def __init__(self,d):
        self.__dict__.update(d)
 
@cython.final
cdef class FlowPathCollection(list):
    
    cpdef update_existence(self, Core):
        """
        A function to update the pointers for the flow path as well as check existence
        of the states on either end of the path
        
        This is required whenever the existence of any of the CV or tubes 
        changes.  Calling this function will update the pointers to the states
        
        Parameters
        ----------
        Core : PDSimCore instance or derived class thereof
        
        """
        cdef FlowPath FP
        cdef dict Tube_Nodes = Core.Tubes.get_Nodes()
        cdef list exists_keys = Core.CVs.exists_keys
        cdef list flow_paths = self
        
        # Set some things in the class that are invariant to save time 
        self.Nexists = len(exists_keys)
        self.N = len(self)
        self.omega = Core.omega
        
        for FP in flow_paths:
            ## Update the pointers to the states for the ends of the flow path
            if FP.key1 in exists_keys:
                CV = Core.CVs[FP.key1]
                FP.State1=CV.State
                FP.key1_exists = True
                FP.ikey1 = exists_keys.index(FP.key1)
            elif FP.key1 in Tube_Nodes:
                FP.State1 = Tube_Nodes[FP.key1]
                FP.key1_exists = False
                if Core.Tubes[FP.key1].key1 == FP.key1:
                    FP.ikey1 = Core.Tubes[FP.key1].i1
                else:
                    FP.ikey1 = Core.Tubes[FP.key1].i2
            else:
                FP.exists=False
                FP.key1_exists = False
                FP.mdot = 0.0
                #Doesn't exist, go to next flow
                continue
            
            if FP.key2 in exists_keys:
                CV = Core.CVs[FP.key2]
                FP.State2=CV.State
                FP.key2_exists = True
                FP.ikey2 = exists_keys.index(FP.key2)
            elif FP.key2 in Tube_Nodes:
                FP.State2=Tube_Nodes[FP.key2]
                FP.key2_exists = False
                if Core.Tubes[FP.key2].key1 == FP.key2:
                    FP.ikey2 = Core.Tubes[FP.key2].i1
                else:
                    FP.ikey2 = Core.Tubes[FP.key2].i2
            else:
                FP.exists=False
                FP.key2_exists = False
                FP.mdot = 0.0
                #Doesn't exist, go to next flow
                continue
            
            #Made it this far, so both states exist
            FP.exists = True
    
    cpdef calculate(self, arraym harray, arraym parray, arraym Tarray):
        """
        Run the code for each flow path to calculate the flow rates
        
        Parameters
        ----------
        harray : :class:`arraym <PDSim.misc.datatypes.arraym>' instance
            arraym that maps index to enthalpy - CVs+Tubes
        parray : :class:`arraym <PDSim.misc.datatypes.arraym>' instance
            arraym that maps index to pressure - CVs+Tubes
        Tarray : :class:`arraym <PDSim.misc.datatypes.arraym>' instance
            arraym that maps index to temperature - CVs+Tubes
        """
        cdef FlowPath FP
        cdef int i
                
        for i in range(self.N):
            FP = self[i]
            if FP.exists:
                FP.calculate(harray, parray, Tarray)
            else:
                FP.edot = 0.0
        
    @cython.cdivision(True)
    cpdef sumterms(self, arraym summerdT, arraym summerdm):
        """
        Sum all the mass flow and mdot*h for each CV in existence at a given 
        step for the derivatives in the ODE solver
        
        Meant to be called by PDSimCore.derivs()
        
        Parameters
        ----------
        summerdT : :class:`arraym <PDSim.misc.datatypes.arraym>' instance
        summerdm : :class:`arraym <PDSim.misc.datatypes.arraym>' instance
            
        """
        cdef double mdot, h_up
        cdef int I_up,I_down
        cdef FlowPath Flow
        cdef int i
        
        summerdm.fill(0.0)
        summerdT.fill(0.0)
        
        #Loop over the flowpaths
        for i in range(self.N):
            
            Flow = self[i]
            
            #One of the chambers doesn't exist if it doesn't have a mass flow term
            if not Flow.exists:
                continue
            else:
                #Do these once for each flow path to cut down on lookups
                mdot = Flow.mdot
                h_up = Flow.h_up
            
            #Flow must exist then
            
            #If the upstream node is a control volume 
            if Flow.key_up_exists:
                #Flow is leaving the upstream control volume
                summerdm.data[Flow.ikey_up] -= mdot/self.omega
                summerdT.data[Flow.ikey_up] -= mdot/self.omega*h_up
                
            #If the downstream node is a control volume
            if Flow.key_down_exists:
                #Flow is entering the downstream control volume
                summerdm.data[Flow.ikey_down] += mdot/self.omega
                summerdT.data[Flow.ikey_down] += mdot/self.omega*h_up
    
    cpdef get_deepcopy(self):
        """
        Using this method, the link to the mass flow function is broken
        """
        return [Flow.get_deepcopy() for Flow in self]

cdef class FlowPath(object):
    
    def __init__(self, key1='', key2='', MdotFcn=None, MdotFcn_kwargs={}):
        """
        
        Parameters
        ----------
        key1 : string
            The key for the first flow node connected to this path
        key2 : string
            The key for the second flow node connected to this path
        MdotFcn : function
            Two options, either an instance of :class:`FlowFunction <PDSim.flow.FlowFunction>', or a 
            function with a prototype like ``f(double A,FlowPath FP, **kwargs)``.  
            See also :class:`FlowFunction <PDSim.flow.FlowFunction>'.  
            
            Any function
            passed in for ``MdotFcn`` will be wrapped into an instance of 
            :class:`FlowFunction <PDSim.flow.FlowFunction>'.  Using an instance of
            :class:`FlowFunction <PDSim.flow.FlowFunction>' is more computationally 
            efficient because the Cython code doesn't need to pass back through 
            the python level and can all stay at the C/C++ level.
            
        MdotFcn_kwargs : dictionary
            A dictionary of terms that will be passed along to the call to 
            ``MdotFcn`` when it is called
        
        """
        self.key1 = key1
        self.key2 = key2
        
        # You are passing in a pre-wrapped function - will be nice and fast since
        # all the calls will stay at the C++ layer
        if isinstance(MdotFcn, FlowFunction):
            self.MdotFcn = MdotFcn
        else:
            # Add the bound method in a wrapper - this will keep the calls at
            # the Python level which will make them nice to deal with but slow
            self.MdotFcn = PyFlowFunctionWrapper(MdotFcn, MdotFcn_kwargs)
        
        self.MdotFcn_str = str(MdotFcn)
            
    cpdef dict __cdict__(self, AddStates=False):
        """
        Returns a dictionary with all the terms that are defined at the Cython level
        """
        
        cdef list items
        cdef list States
        cdef dict d={}
        cdef bytes item
        
        items=['mdot','h_up','h_down','T_up','p_up','p_down','key_up','key_down','key1','key2','Gas','exists']
        for item in items:
            d[item]=getattr(self,item)
        
        if AddStates:
            States=[('State1',self.State1),('State2',self.State2),('State_up',self.State_up),('State_down',self.State_down)]
            for k,State in States:
                if State is not None:
                    d[k]=State.copy()
        return d
    
    cpdef FlowPath get_deepcopy(self):
        cdef FlowPath FP = FlowPath.__new__(FlowPath)
        FP.Gas = self.Gas
        FP.exists = self.exists
        FP.mdot = self.mdot
        FP.edot = self.edot
        FP.h_up = self.h_up
        FP.T_up = self.T_up
        FP.p_up = self.p_up
        FP.p_down = self.p_down
        FP.key_up = self.key_up
        FP.key_down = self.key_down
        FP.key_up_exists = self.key_up_exists
        FP.key_down_exists = self.key_down_exists
        FP.A = self.A
        FP.key1_exists = self.key1_exists
        FP.key2_exists = self.key2_exists
#         if self.exists:
#             FP.State_up = self.State_up.copy()
#             FP.State_down = self.State_down.copy()
        return FP
        
    cpdef calculate(self, arraym harray, arraym parray, arraym Tarray):
        """
        Calculate all of the flow paths
        
        Parameters
        ----------
        harray : :class:`arraym <PDSim.misc.datatypes.arraym>' instance of enthalpies of CV+Tubes
        parray : :class:`arraym <PDSim.misc.datatypes.arraym>' instance of pressures of CV+Tubes
        Tarray : :class:`arraym <PDSim.misc.datatypes.arraym>' instance of temperatures of CV+Tubes
        """
        cdef FlowFunction FF
        cdef double p1 = parray.data[self.ikey1], p2 = parray.data[self.ikey2]
        
        if p1 > p2:
            # The pressure in chamber 1 is higher than chamber 2
            # and thus the flow is from chamber 1 to 2
            self.key_up = self.key1
            self.key_down = self.key2
            self.State_up = self.State1
            self.State_down = self.State2
            self.T_up = Tarray.data[self.ikey1]
            self.h_up = harray.data[self.ikey1]
            self.h_down = harray.data[self.ikey2]
            self.p_up = p1
            self.p_down = p2
            self.key_up_exists = self.key1_exists
            self.key_down_exists = self.key2_exists
            self.key_up_Index = self.key1Index
            self.key_down_Index = self.key2Index
            self.ikey_up = self.ikey1
            self.ikey_down = self.ikey2
            self.Gas = self.State1.Fluid
        else:
            self.key_up = self.key2
            self.key_down = self.key1
            self.State_up = self.State2
            self.State_down = self.State1
            self.T_up = Tarray.data[self.ikey2]
            self.h_up = harray.data[self.ikey2]
            self.h_down = harray.data[self.ikey1]                  
            self.p_up = p2
            self.p_down = p1
            self.key_up_exists = self.key2_exists
            self.key_down_exists = self.key1_exists
            self.key_up_Index = self.key2Index
            self.key_down_Index = self.key1Index
            self.ikey_up = self.ikey2
            self.ikey_down = self.ikey1
            self.Gas = self.State2.Fluid
            
        FF = self.MdotFcn
        self.mdot = FF.call(self)
        
        self.edot = abs(self.mdot*((self.h_up - self.h_down)-298.15*(self.State_up.get_s()-self.State_down.get_s())  ))
        
    def __reduce__(self):
        return rebuildFlowPath,(self.__getstate__(),)
        
    def __getstate__(self):
        d={}
        d['MdotFcn']=cPickle.dumps(self.MdotFcn)
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
        
def rebuildFlowPath(d):
    FP = FlowPath()
    FP.MdotFcn = cPickle.loads(d.pop('MdotFcn'))
    for item in d:
        setattr(FP,item,d[item])
    return FP
    