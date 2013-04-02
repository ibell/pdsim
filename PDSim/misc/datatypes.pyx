import cython
cimport cython

import numpy as np
cimport numpy as np
    
cimport cpython.array
from libc.stdlib cimport calloc, free, realloc
from libc.string cimport memcpy
from cpython cimport bool

cimport cython

cdef class AnnotatedValue(object):
    def __init__(self, str key, object value, str annotation, str units):
        self.key = key
        self.value = value
        self.annotation = annotation
        self.units = units

cdef class Collector(object):
    """
    The collector class is a thin wrapper around a list with compact semantics for adding values::
    
        C  = Collector()
        C << 1
        C << 3
        
    And the values can be obtained as a numpy array by::
    
        C.v()
        
    Or as a list from the attribute:: 
    
        C.vec
        
    """
    def __init__(self):
        self.vec = []
        
    def __lshift__(self, other):
        self.vec.append(other)
        
    def __getitem__(self, item):
        return self.vec[item]
        
    def __repr__(self):
        return str(self.vec)
    
    cpdef v(self, int ndmin = 2):
        """
        Get a numpy array of the Collector values
        
        Internally, the call looks like::
        
            return np.array(self.vec, ndmin = ndmin)
            
        Parameters
        ----------
        ndmin, int
            Minimum number of axes
        """
        return np.array(self.vec, ndmin = ndmin)
    
@cython.final
cdef class arraym(object):
    """
    A thin wrapper around a low-level c-array with pythonic semantics
    
    Implements the element-wise operators +,-,/,* between arraym instances and other iterables and constants
    
    Note: The only divide operator implemented is the truediv operator of python 3.x .  To get the ``/`` operator to work in python 2.x, you need to do ``from __future__ import division`` in your header
    """
    
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
            if isinstance(data,(float,int)):
                #It is a single value, wrap it in a list
                self.N = 1
                data = [data]
            else:
                #Its an item with a length
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
            else:
                #Now it must either be an iterable or a failure
                try:
                    for i,el in enumerate(data):
                        self.data[i] = el
                except TypeError:
                    raise TypeError("Sorry but you provided a type to arraym that doesn't work.  Good types are arraym, numpy arrays, or any iterable.")
        else:
            self.data = NULL

    def __cinit__(self):
        self.N = 0
        self.data = NULL
        
    cpdef set_size(self, int N):
        """
        Set the size of the internal array, initialized to zeros
        """
        #If a zero or negative length passed in, don't allocate memory, but set length flag 
        if N <= 0:
            self.N = 0
            return
        
        if self.data == NULL:
            #Allocate the memory for the array that will be used internally
            self.data = <double *> calloc(N, sizeof(double))
            self.N = N
            
    cdef void set_data(self, double *data, int N):
        """
        Parameters
        ----------
        data : double *
            The data to be set in this array
        N : int
            The size of the passed data
        """
        if self.data == NULL:
            #Allocate the memory for the array that will be used internally
            self.data = <double *> calloc(N, sizeof(double))
            self.N = N
        elif not self.N == N:
            raise ValueError('Memory already allocated for arraym, but sizes of arraym ('+str(self.N)+') and data ('+str(N)+') do not match')
        memcpy(self.data,data,N*sizeof(double))
    
    cpdef dealloc(self):
        """ Clean up the memory we allocated. Not generally needed to be called, python will do the cleanup automatically """
        if not self.data == NULL:
            free(self.data)
            self.data = NULL
            self.N = 0
    
    def __dealloc__(self):
        #Clean up the memory we allocated
        if not self.data == NULL:
            free(self.data)
            self.data = NULL
            self.N = 0
          
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
            try:
                #Try to make an iterator out of y
                iterator = iter(y)
            except TypeError:
                # not iterable - int, float, etc.
                
                # Cast y to a double
                yd = (<double>y)
                # Add on the other array values
                for i in range(N):
                    zdata[i] += yd
            else:
                # iterable - list, tuple, numpy array, etc.
                for i in range(N):
                    zdata[i] += y[i]
        
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
            try:
                #Try to make an iterator out of y
                iterator = iter(y)
            except TypeError:
                # not iterable - int, float, etc.
                
                # Cast to a double
                yd = (<double>y)
                # Multiply by the other value
                for i in range(N):
                    zdata[i] *= yd
            else:
                # iterable - list, tuple, numpy array, etc.
                #No type introspection possible
                for i in range(N):
                    zdata[i] *= y[i]
            
        return z
    
    def __truediv__(x, y):
        cdef int i, N
        cdef bint isarray_x, isarray_y
        cdef double *zdata, *ydata
        cdef double yd,xd
        cdef arraym z
        
        isarray_x = isinstance(x, arraym)
        isarray_y = isinstance(y, arraym)

        #Both x and y are arraym instances
        if isarray_x & isarray_y:
            check_dims(x, y)
            N = (<arraym>x).N
            z = (<arraym>x).copy()
            zdata = (<arraym>z).data
            ydata = (<arraym>y).data
            # Add on the other array values
            for i in range(N):
                zdata[i] /= ydata[i]
                
        #One of x and y is an arraym
        elif isarray_x != isarray_y:
            if isarray_y:
                N = (<arraym>y).N
                z = (<arraym>y).copy()
                zdata = (<arraym>z).data
                if isinstance(x,(int,float)):
                    # Cast lhs to a double and rhs to a double*
                    xd = (<double>x)
                    ydata = (<arraym>y).data
                    # Add on the other array values
                    for i in range(N):
                        zdata[i] = xd/ydata[i]
                else:
                    #Hopefully it is an iterable
                    for i in range(len(x)):
                        z[i] = x[i]/y[i]                    
            else:
                N = (<arraym>x).N
                z = (<arraym>x).copy()
                zdata = (<arraym>z).data
                if isinstance(y,(int,float)):
                    # Cast rhs to a double
                    yd = <double> y
                    # Add on the other array values
                    for i in range(N):
                        zdata[i] /= yd
                else:
                    #Hopefully it is an iterable
                    for i in range(len(x)):
                        z[i] = x[i]/y[i]
                    
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
                
        #One of x and y is an arraym
        elif isarray_x != isarray_y:
            if isarray_y:
                N = (<arraym>y).N
                z = (<arraym>y).copy()
                zdata = (<arraym>z).data
                if isinstance(x,(int,float)):
                    # Cast lhs to a double and rhs to a double*
                    xd = (<double>x)
                    ydata = (<arraym>y).data
                    # Add on the other array values
                    for i in range(N):
                        zdata[i] = xd - ydata[i]
                else:
                    #Hopefully it is an iterable
                    for i in range(len(x)):
                        z[i] = x[i] - y[i]                    
            else:
                N = (<arraym>x).N
                z = (<arraym>x).copy()
                zdata = (<arraym>z).data
                if isinstance(y,(int,float)):
                    # Cast rhs to a double
                    yd = <double> y
                    # Add on the other array values
                    for i in range(N):
                        zdata[i] -= yd
                else:
                    #Hopefully it is an iterable
                    for i in range(len(x)):
                        z[i] = x[i] - y[i]
            
        return z
    
    cpdef arraym copy(self):
        """
        Return a copy of the instance
        """
        cdef arraym arr = arraym.__new__(arraym)
        arr.set_data(self.data, self.N)
        return arr
    
    def __setitem__(self,int i, double y):
        self.data[i]=y
        
    @cython.returns(double)
    def __getitem__(self, int i):
        return self.data[i]
    
    cpdef double get_index(self, int i) except *:
        """
        Get the value at the given index
        
        Parameters
        ----------
        i : int
            The index in the array
        """
        if i < self.N:
            return self.data[i]
        else:
            raise ValueError('Your index [{i:d}] is out of range [{N:d}]'.format(i=i,N=self.N))
    
    cpdef double set_index(self, int i, double val) except *:
        """
        Set the value at the given index
        
        Parameters
        ----------
        i : int
            The index in the array
        """
        if i < self.N:
            self.data[i] = val
        else:
            raise ValueError('Your index [{i:d}] is out of range [{N:d}]'.format(i=i,N=self.N))
        
    cpdef fill(self, double fillval):
        """ Fill the array with the given value """
        for i in range(self.N):
            self.data[i] = fillval
    
    cdef arraym slice(self, int i, int j):
        cdef int k
        cdef arraym arr = arraym()
        if j < i:
            raise IndexError('Indices must be increasing')
        if j == i:
            raise IndexError('Length of slice must be greater than 1')
        if j > self.N:
            raise IndexError('End of slice out of bounds. Length of arraym is '+str(self.N)+' and requested end is '+str(j))
        
        arr.set_size(j-i)
        memcpy(arr.data,self.data+i,(j-i)*sizeof(double))
        return arr
    
    cpdef extend(self, arraym array2):
        """
        Extend this arraym instance with another arraym instance
        
        Parameters
        ----------
        array2 : :class:`arraym <PDSim.misc.datatypes.arraym>` instance
            The arraym to be appended to this arraym instance 
        """
        cdef double* new_data
        cdef int N = array2.N + self.N
        
        #Only extend if there is something in the extension array
        if N > self.N:
            #Reallocate the array to extend its length
            new_data = <double*>realloc(self.data, N*sizeof(double))
            #Copy into the new array
            memcpy(new_data+self.N, array2.data, array2.N*sizeof(double))
            #Free the old array
            free(self.data)
            #Make self.data point to the newly allocated array
            self.data = new_data
            #Set the length
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

cpdef arraym empty_arraym(int N):
    """
    A convenience function to return an arraym with the given size initialized to zero
    
    Parameters
    ----------
    N: int
        Size of the arraym to return initialized to zero
    """
    cdef arraym arr = arraym()
    arr.set_size(N)
    return arr


###############################################################################
###############################################################################
##          LISTM - Enhanced list with element-wise operators                ##
###############################################################################
###############################################################################


cdef class listm(list):
    """
    See http://docs.cython.org/src/userguide/special_methods.html
    """
    def __add__(self,y):
        cdef int i,N
        cdef bool isarray_x,isarray_y
        
        isarray_x = isinstance(self,listm)
        isarray_y = isinstance(y,listm)
        
        if isinstance(self,listm):
            N=len(self)
            if isinstance(y,int) or isinstance(y,float):
                return listm([self[i]*y for i in range(N)])
            else:
                return listm([self[i]*y[i] for i in range(N)])
        else:
            ### it is backwards, self is something else, y is a listm
            N=len(y)
            if isinstance(self,(int,float)) and not isinstance(y,(int,float)):
                self,y = y,self
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
    
    def copy(self):
        """
        Return a copy of the listm instance
        """
        return listm(self[:])
          
def rebuildListm(d):
    return listm(d['data'])