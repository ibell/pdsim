import cython
cimport cython
    
cimport cpython.array
from libc.stdlib cimport calloc, free
from cpython cimport bool
#cdef class arraym(object):
#    def __init__(self, data = None):
#        """
#        Parameters
#        ----------
#        If data not included, make an empty arraym
#        """
#        cdef int N = len(data)
#        cdef bytes format = bytes("d")
#        cdef tuple shape = (N,)
#        cdef long item_size = sizeof(double)
#        if data is not None:
#            
#            #Allocate the memory for the array that will be used internally
#            self.data = cvarray(shape, item_size, format)
#            
#            #Low-level copy of the data
#            self.data[:] = data
#        
#    cpdef copy(self):   
#        cdef double* ptr = <double*> malloc(sizeof(double) * 100)
#        free(ptr)
##        cdef double data[100]
##        for i in range(100):
##            data[i]=4
#        #cdef double[:] vdata = data
#        #vdata[:] = 2
#        pass
#        #cdef cvarray data = cvarray(shape=(100,), itemsize=sizeof(double), format="d")
#        #Allocate the memory for the array that will be used internally
#        #cdef cvarray data = self.data[:]#cvarray(shape=(self.data.shape[0],), itemsize=sizeof(double), format="d")
#        
#         
#    def __add__(x, y):
#        cdef int i, N, isarray_x, isarray_y
#        cdef double[:] xdata
#        cdef double[:] ydata
#        cdef double[:] zdata
#        isarray_x = isinstance(x, arraym)
#        isarray_y = isinstance(y, arraym)
#
#        if isarray_x & isarray_y:
#            xdata = (<arraym>x).vdata
#            ydata = (<arraym>y).vdata
#            zdata = xdata.copy()
#            for i in range(xdata.shape[0]):
#                zdata[i] += ydata[i]
#            return arraym(zdata)
#    
#    def __iadd__(self,y):
#        cdef int i, isarray_y
#        cdef double[:] ydata
#        isarray_y = isinstance(y, arraym)
#
#        if isarray_y:
#            ydata = (<arraym>y).vdata
#            for i in range(self.vdata.shape[0]):
#                self.vdata[i] += ydata[i]
#        
#    def __repr__(self):
#        return str(list(self.data))

#cdef class arraym(cvarray):
#    def __init__(self, data):
#        cdef int i
#        carray = cvarray(shape = (len(data),), itemsize=sizeof(double), format='d')
#        for i in range(len(data)):
#            carray[i]=data[i]
#        print carray[0]

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