import copy
import cPickle
from flow_models import PyFlowFunctionWrapper
from PDSim.misc._listmath import listm

cpdef sumterms_given_CV(bytes key, list Flows):
    cdef FlowPath Flow
    cdef double summer_mdot,summer_mdoth
    
    summer_mdot = 0.0
    summer_mdoth = 0.0
    
    for Flow in Flows:

        if abs(Flow.mdot)<1e-12:
            continue
    
        if Flow.key_down==key:
            summer_mdot+=Flow.mdot
            summer_mdoth+=Flow.mdot*Flow.h_up
        
        elif Flow.key_up==key:
            summer_mdot-=Flow.mdot
            summer_mdoth-=Flow.mdot*Flow.h_up
            
        else:
            continue
        
    return summer_mdot, summer_mdoth

class struct(object):
    def __init__(self,d):
        self.__dict__.update(d)
 
cdef class FlowPathCollection(list):
    def __init__(self,d=None):
        if d is not None and isinstance(d,dict):
            self.__dict__.update(d)
    
    cpdef calculate(self, Core, dict hdict):
        """
        Run the code for each flow path to calculate the flow rates & Update the states for the CVs and the Tubes.
        
        This is required whenever the existence of any of the CV or tubes 
        changes.  Calling this function will update the pointers to the states
        
        Parameters
        ----------
        Core : PDSimCore instance or derived class thereof
        hdict : dictionary
            Maps CV key to enthalpy
        """
        cdef FlowPath FP
        cdef dict Tubes_Nodes
        cdef list exists_keys
        exists_keys = Core.CVs.exists_keys
        Tubes_Nodes = Core.Tubes.Nodes
                
        for FP in self:
            ## Update the pointers to the states for the ends of the flow path
            if FP.key1 in exists_keys:
                FP.State1=Core.CVs[FP.key1].State
            elif FP.key1 in Tubes_Nodes:
                FP.State1=Tubes_Nodes[FP.key1]
            else:
                FP.mdot=0.0
                #Doesn't exist, go to next flow
                continue                   
            
            if FP.key2 in exists_keys:
                FP.State2=Core.CVs[FP.key2].State
            elif FP.key2 in Tubes_Nodes:
                FP.State2=Tubes_Nodes[FP.key2]
            else:
                FP.mdot=0.0
                #Doesn't exist, go to next flow
                continue
            
            FP.calculate(hdict)
        
    @property
    def N(self):
        return self.__len__()
        
    cpdef sumterms_helper(self, list exists_keys, double omega):
        cdef list summerdm, summerdT
        cdef double mdot, h_up
        cdef int I_up,I_down
        cdef bytes key_up, key_down
        cdef FlowPath Flow
        
        summerdm = [0.0 for _dummy in range(len(exists_keys))]
        summerdT = summerdm[:]
        
        for Flow in self:
            #One of the chambers doesn't exist if it doesn't have a mass flow term
            if Flow.mdot==0:
                continue
            else:
                mdot=Flow.mdot
                h_up = Flow.h_up
            
            #Flow must exist then
            key_up=Flow.key_up    
            if key_up in exists_keys:
                #The upstream node is a control volume
                I_up=exists_keys.index(key_up)
                #Flow is leaving the upstream control volume
                summerdm[I_up]-=mdot/omega
                summerdT[I_up]-=mdot/omega*h_up
                
            key_down=Flow.key_down
            if key_down in exists_keys:
                #The downstream node is a control volume                
                I_down=exists_keys.index(key_down)
                #Flow is entering the downstream control volume
                summerdm[I_down]+=mdot/omega
                summerdT[I_down]+=mdot/omega*h_up
    
        return summerdm, summerdT

    
    
    cpdef sumterms(self,Core):
        
        #Call the Cython-version of the summer
        summerdm,summerdT = self.sumterms_helper(Core.CVs.exists_keys, Core.omega)
        
        return listm(summerdT),listm(summerdm)
    
    def __reduce__(self):
        return rebuildFPC,(self[:],)
    
    cpdef deepcopy(self):
        newFPC = FlowPathCollection()
        newFPC.__dict__.update(copy.deepcopy(self.__dict__))
        return newFPC
    
    cpdef get_deepcopy(self):
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

def rebuildFPC(d):
    """
    Used with cPickle to recreate the (empty) flow path collection class
    """
    FPC = FlowPathCollection()
    for Flow in d:
        FPC.append(Flow)
    return FPC

cdef class FlowPath(object):
    
    def __init__(self, key1='', key2='', MdotFcn=None, MdotFcn_kwargs={}):
        self.key1 = key1
        self.key2 = key2
        
        # You are passing in a pre-wrapped function - will be nice and fast since
        # all the calls will stay at the C++ layer
        if isinstance(MdotFcn, FlowFunctionWrapper):
            self.MdotFcn = MdotFcn
        else:
            # Add the bound method in a wrapper - this will keep the calls at
            # the Python level which will make them nice to deal with but slow
            self.MdotFcn = PyFlowFunctionWrapper(MdotFcn, MdotFcn_kwargs)
            
    cpdef dict __cdict__(self, AddStates=False):
        """
        Returns a dictionary with all the terms that are defined at the Cython level
        """
        cdef list items=['mdot','h_up','T_up','p_up','p_down','key_up','key_down'
                         ,'key1','key2','Gas']
        cdef list States=[('State1',self.State1),('State2',self.State2),
                          ('State_up',self.State_up),('State_down',self.State_down)]
        cdef dict dic={}
        for item in items:
            dic[item]=getattr(self,item)
        if AddStates==True:
            for k,State in States:
                if State is not None:
                    dic[k]=State.copy()
        return dic
    
    cpdef FlowPath get_deepcopy(self):
        cdef dict d = self.__cdict__()
        cdef FlowPath FP = FlowPath()
        
        for k,v in d.iteritems():
            setattr(FP,k,v)
        return FP
        
    cpdef calculate(self, dict hdict = None):
        """
        calculate
        """
        cdef double p1,p2
        cdef FlowFunctionWrapper FW
        #The states of the chambers
        p1=self.State1.get_p()
        p2=self.State2.get_p()
        
        if p1 > p2:
            # The pressure in chamber 1 is higher than chamber 2
            # and thus the flow is from chamber 1 to 2
            self.key_up=self.key1
            self.key_down=self.key2
            self.State_up=self.State1
            self.State_down=self.State2
            self.T_up=self.State1.get_T()
            if hdict is None:
                self.h_up=self.State1.get_h()
            else:
                self.h_up = hdict[self.key_up]
            self.p_up=p1
            self.p_down=p2
            
            self.Gas=self.State1.Fluid
        else:
            self.key_up=self.key2
            self.key_down=self.key1
            self.State_up=self.State2
            self.State_down=self.State1
            self.T_up=self.State2.get_T()
            if hdict is None:
                self.h_up=self.State2.get_h()
            else:
                self.h_up = hdict[self.key_up]
            self.p_up=p2
            self.p_down=p1
            
            self.Gas=self.State2.Fluid
            
        FW = self.MdotFcn
        self.mdot = FW.call(self)
        
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
    