import cython
cimport cython

import numpy as np
cimport numpy as np
    
cimport cpython.array
from libc.stdlib cimport calloc, free
from libc.string cimport memcpy
from cpython cimport bool

cdef class arraym(object):
    def __init__(self, data = None):
        """
        data : list, array.array, numpy array, etc.
        
        Notes
        -----
        If a numpy array is provided, the numpy buffer is used internally to access the data
        Otherwise, as long as the iterable contains floating point values it should work
        """
        cdef int i
        cdef double el
        cdef np.ndarray[np.float_t, ndim = 1] npdata
        
        if data is not None:
            self.N = len(data)
            #Allocate the memory for the array that will be used internally
            self.data = <double *> calloc(self.N, sizeof(double))
            
            #If a numpy array use the buffering interface
            if isinstance(data, np.ndarray):
                npdata = data
                for i in range(self.N):
                    self.data[i] = npdata[i]
            elif isinstance(data, list):
                for i,el in enumerate(data):
                    self.data[i] = el
        else:
            self.data = NULL
            
    cdef set_data(self, double *data, int N):
        cdef int i
        if self.data is NULL:
            #Allocate the memory for the array that will be used internally
            self.data = <double *> calloc(N, sizeof(double))
            self.N = N
        memcpy(self.data,data,N*sizeof(double))
            
    def __dealloc__(self):
        #Clean up the memory we allocated
        if self.data is not NULL:
            free(self.data)
          
    def __add__(x, y):
        cdef int i, N
        cdef bint isarray_x, isarray_y
        cdef double *zdata,*ydata,*xdata
        cdef double yd
        cdef arraym z
        
        isarray_x = isinstance(x, arraym)
        isarray_y = isinstance(y, arraym)

        if isarray_x & isarray_y:
            N = (<arraym>x).N
            z = (<arraym>x).copy()
            zdata = (<arraym>z).data
            ydata = (<arraym>y).data
            # Add on the other array values
            for i in range(N):
                zdata[i] += ydata[i]
        elif isarray_x != isarray_y:
            if isarray_y:
                x,y = y,x
            N = (<arraym>x).N
            z = (<arraym>x).copy()
            zdata = (<arraym>z).data
            # Cast to a double
            yd = <double> y
            # Add on the other array values
            for i in range(N):
                zdata[i] += yd
        
        return z
    
    def __mul__(x, y):
        cdef int i, N
        cdef bint isarray_x, isarray_y
        cdef double *zdata,*ydata
        cdef double yd
        cdef arraym z
        
        isarray_x = isinstance(x, arraym)
        isarray_y = isinstance(y, arraym)

        if isarray_x & isarray_y:
            N = (<arraym>x).N
            z = (<arraym>x).copy()
            zdata = (<arraym>z).data
            ydata = (<arraym>y).data
            for i in range(N):
                zdata[i] *= ydata[i]
        elif isarray_x != isarray_y:
            if isarray_y:
                x,y = y,x
            N = (<arraym>x).N
            z = (<arraym>x).copy()
            zdata = (<arraym>z).data
            # Cast to a double
            yd = <double> y
            # Add on the other array values
            for i in range(N):
                zdata[i] *= yd
            
        return z
    
    def __truediv__(x, y):
        cdef int i, N
        cdef bint isarray_x, isarray_y
        cdef double *zdata, *ydata
        cdef double yd,xd
        cdef arraym z
        
        isarray_x = isinstance(x, arraym)
        isarray_y = isinstance(y, arraym)

        if isarray_x & isarray_y:
            N = (<arraym>x).N
            z = (<arraym>x).copy()
            zdata = (<arraym>z).data
            ydata = (<arraym>y).data
            # Add on the other array values
            for i in range(N):
                zdata[i] /= ydata[i]
        elif isarray_x != isarray_y:
            if isarray_y:
                N = (<arraym>y).N
                z = (<arraym>y).copy()
                zdata = (<arraym>z).data
                 # Cast lhs to a double and rhs to a double*
                xd = <double> x
                ydata = (<arraym>y).data
                # Add on the other array values
                for i in range(N):
                    zdata[i] = xd/ydata[i]
            else:
                N = (<arraym>x).N
                z = (<arraym>x).copy()
                zdata = (<arraym>z).data
                 # Cast rhs to a double
                yd = <double> y
                # Add on the other array values
                for i in range(N):
                    zdata[i] /= yd
            
        return z
    
    def __sub__(x, y):
        cdef int i, N
        cdef bint isarray_x, isarray_y
        cdef double *zdata, *ydata
        cdef double yd,xd
        cdef arraym z
        
        isarray_x = isinstance(x, arraym)
        isarray_y = isinstance(y, arraym)

        if isarray_x & isarray_y:
            N = (<arraym>x).N
            z = (<arraym>x).copy()
            zdata = (<arraym>z).data
            ydata = (<arraym>y).data
            # Add on the other array values
            for i in range(N):
                zdata[i] -= ydata[i]
        elif isarray_x != isarray_y:
            if isarray_y:
                N = (<arraym>y).N
                z = (<arraym>y).copy()
                zdata = (<arraym>z).data
                 # Cast lhs to a double and rhs to a double*
                xd = <double> x
                ydata = (<arraym>y).data
                # Add on the other array values
                for i in range(N):
                    zdata[i] = xd - ydata[i]
            else:
                N = (<arraym>x).N
                z = (<arraym>x).copy()
                zdata = (<arraym>z).data
                 # Cast rhs to a double
                yd = <double> y
                # Add on the other array values
                for i in range(N):
                    zdata[i] -= yd
            
        return z
        
    cdef arraym copy(self):
        cdef arraym arr = arraym()
        arr.set_data(self.data, self.N)
        return arr
    
    def __iter__(self):
        for i in range(self.N):
            yield float(self.data[i])
        
    def __repr__(self):
        return str(list(self))

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