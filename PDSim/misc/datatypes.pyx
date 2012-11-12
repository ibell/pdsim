import cython
cimport cython

import numpy as np
cimport numpy as np
    
cimport cpython.array
from libc.stdlib cimport calloc, free, realloc
from libc.string cimport memcpy
from cpython cimport bool

cimport cython

@cython.final
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
            
            #If it is already an arraym instance, do a low-level copy of the data
            if isinstance(data,arraym):
                memcpy(self.data,(<arraym>data).data,self.N*sizeof(double))
            #If a numpy array use the buffering interface
            elif isinstance(data, np.ndarray): 
                npdata = data
                for i in range(self.N):
                    self.data[i] = npdata[i]
            #If it is an array.array, use the buffer interface
            #elif isinstance(data, array.array):
            #    vdata = data
            #    self.data[:] = vdata
            elif isinstance(data, list):
                for i,el in enumerate(data):
                    self.data[i] = el
        else:
            self.data = NULL
            
    cpdef set_size(self, int N):
        """
        Set the size of the internal array, initialized to zeros
        """
        if self.data is NULL:
            #Allocate the memory for the array that will be used internally
            self.data = <double *> calloc(N, sizeof(double))
            self.N = N
            
    cdef void set_data(self, double *data, int N):
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
        cdef double *zdata,*ydata
        cdef double yd
        cdef arraym z
        
        isarray_x = isinstance(x, arraym)
        isarray_y = isinstance(y, arraym)

        if isarray_x & isarray_y:
            check_dims(x, y)
            N = (<arraym>x).N
            z = (<arraym>x).copy()
            zdata = z.data
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
            yd = (<double>y)
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
            check_dims(x, y)
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
            yd = (<double>y)
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
            check_dims(x, y)
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
                xd = (<double>x)
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
            check_dims(x, y)
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
                xd = (<double>x)
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
        
    cpdef arraym copy(self):
        cdef arraym arr = arraym()
        arr.set_data(self.data, self.N)
        return arr
    
    def __setitem__(self,int i, double y):
        self.data[i]=y
        
    @cython.returns(double)
    def __getitem__(self, int i):
        return self.data[i]
    
    cdef arraym slice(self, int i, int j):
        cdef int k
        cdef arraym arr = arraym()
        if j < i:
            raise IndexError('Indices must be increasing')
        if j == i:
            raise IndexError('Length of slice must be greater than 1')
        if j > self.N:
            raise IndexError('End of slice out of bounds. Length of arraym is '+str(self.N))
        
        arr.set_size(j-i)
        memcpy(arr.data,self.data+i,(j-i)*sizeof(double))
        return arr
    
    cpdef extend(self, arraym array2):
        cdef int N = array2.N + self.N
        self.data = <double*>realloc(self.data, N*sizeof(double))        
        memcpy(self.data+self.N, array2.data, array2.N*sizeof(double))
        self.N = N
        
    def __getslice__(self, Py_ssize_t i, Py_ssize_t j):
        return self.slice(i,j)
    
    def __iter__(self):
        for i in range(self.N):
            yield float(self.data[i])
        
    def __repr__(self):
        return str(list(self))
        
    def __len__(self):
        return self.N
