import cython
cimport cython

from CoolProp.State import State as StateClass
from CoolProp.State cimport State as StateClass

import numpy as np
cimport numpy as np

from PDSim.flow._flow import _FlowPath
from PDSim.flow._flow cimport _FlowPath

from PDSim.misc._listmath import listm
from PDSim.misc._listmath cimport listm

cpdef list getcol(np.ndarray[np.float_t, ndim=2] mat, int colIndex, list indices):
    cdef int i,j
    cdef list values = []
    for i in indices: 
        values.append(mat[i,colIndex])
    #print 'getcol,i=',i,'indices=',indices,'vals=',values
    return values

cpdef setcol(np.ndarray[np.float_t, ndim=2] mat, int colIndex, list indices, listm x):
    cdef int i,j
    cdef double val
    #print 'setcol,i=',colIndex,'indices=',indices,'vals=',x
    if not len(indices) == len(x):
        print 'indices:',indices, 'x:',x
        raise ValueError('Length of indices [{0:d}] and x [{1:d}] are not the same'.format(len(indices),len(x)))
    for j in range(len(x)):
        i = indices[j]
        val = x[j]
        mat[i,colIndex]=val
    return None

cpdef sumterms_helper(object Flows, list exists_keys, double omega):
    cdef list summerdm, summerdT
    cdef double mdot, h_up
    cdef int I_up,I_down
    cdef bytes key_up, key_down
    cdef _FlowPath Flow
    
    summerdm = [0.0 for _dummy in range(len(exists_keys))]
    summerdT = summerdm[:]
    
    for Flow in Flows:
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
