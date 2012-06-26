import cython
cimport cython
        
cdef class listm(list):
    """
    See http://docs.cython.org/src/userguide/special_methods.html
    """
    def __add__(self,y): 
        cdef int i,N
        
        if isinstance(self,listm):
            N=len(self)
            if isinstance(y,int) or isinstance(y,float):
                return listm([self[i]+y for i in range(N)])
            else:
                return listm([self[i]+y[i] for i in range(N)])
        else:
            ### it is backwards, self is something else, y is a listm
            N=len(y)
            if isinstance(self,int) or isinstance(self,float):
                return listm([y[i]+self for i in range(N)])
            else:
                return listm([self[i]+y[i] for i in range(N)])
                
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
          
def rebuildListm(d):
    return listm(d['data'])