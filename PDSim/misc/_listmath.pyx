import cython
cimport cython

import numpy as np
cimport numpy as np
    
cimport cpython.array
from libc.stdlib cimport calloc, free
from libc.string cimport memcpy
from cpython cimport bool

cdef class listm(list):
    """
    See http://docs.cython.org/src/userguide/special_methods.html
    """
    def __add__(x,y):
        cdef int i,N
        cdef bool isarray_x,isarray_y
        
        isarray_x = isinstance(x,listm)
        isarray_y = isinstance(y,listm)
        
        if isarray_x and isarray_y:
            return listm([x[i]+y[i] for i in range(x.__len__())])
        else:
            ### it is backwards, self is something else, y is a listm
            N=len(y)
            if isinstance(x,(int,float)) and not isinstance(y,(int,float)):
                x,y = y,x
                return listm([y[i]+x for i in range(N)])
            else: 
                return listm([x[i]+y[i] for i in range(N)])
                
    def __mul__(self,y):
        cdef int i,N
        
        if isinstance(self,listm):
            N=len(self)
            if isinstance(y,int) or isinstance(y,float):
                return listm([self[i]*y for i in range(N)])
            else:
                return listm([self[i]*y[i] for i in range(N)])
        else:
            ### it is backwards, self is something else, y is a listm
            N=len(y)
            if isinstance(self,int) or isinstance(self,float):
                return listm([y[i]*self for i in range(N)])
            else:
                return listm([self[i]*y[i] for i in range(N)])
                
    def __truediv__(self,y):
        cdef int i,N
        
        if isinstance(self,listm):
            N=len(self)
            if isinstance(y,int) or isinstance(y,float):
                return listm([self[i]/y for i in range(N)])
            else:
                return listm([self[i]/y[i] for i in range(N)])
        else:
            ### it is backwards, self is something else, y is a listm
            N=len(y)
            if isinstance(self,int) or isinstance(self,float):
                return listm([self/y[i] for i in range(N)])
            else:
                return listm([self[i]/y[i] for i in range(N)])
    
    def __sub__(self,y):
        cdef int i,N
        
        if isinstance(self,listm):
            N=len(self)
            if isinstance(y,int) or isinstance(y,float):
                return listm([self[i]-y for i in range(N)])
            else:
                return listm([self[i]-y[i] for i in range(N)])
        else:
            ### it is backwards, self is something else, y is a listm
            N=len(y)
            if isinstance(self,int) or isinstance(self,float):
                return listm([self-y[i] for i in range(N)])
            else:
                return listm([self[i]-y[i] for i in range(N)])
    
    def __reduce__(self):
        d={}
        d['data']=list(self)
        return rebuildListm,(d,)
    
    def copy(self):
        """
        Return a copy of the listm instance
        """
        return listm(self[:])
          
def rebuildListm(d):
    return listm(d['data'])