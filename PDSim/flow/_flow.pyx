import copy

cpdef tuple sum_flows(bytes key, object Flows):
    cdef _FlowPath Flow
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

cdef class _FlowPathCollection(list):
    pass

cdef class _FlowPath:
    
    cpdef dict __cdict__(self,AddStates=False):
        """
        Returns a dictionary with all the terms that are defined at the Cython level
        """
        cdef list items=['mdot','h_up','T_up','p_up','p_down','key_up','key_down'
                         ,'key1','key2','MdotFcn_kwargs','Gas']
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
    
    cpdef _FlowPath get_deepcopy(self):
        cdef dict d = self.__cdict__()
        cdef _FlowPath FP = _FlowPath()
        
        for k,v in d.iteritems():
            setattr(FP,k,v)
        return FP
        
    cpdef _calculate(self, dict hdict = None):
        """
        calculate
        """
        cdef double p1,p2
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