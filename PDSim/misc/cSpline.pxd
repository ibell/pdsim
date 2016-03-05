from libcpp.vector cimport vector

cdef extern from "spline.h":
        
    cdef cppclass Spline[T,V]:
        
        ## Constructor
        Spline(vector[T] x, vector[V] y) except +ValueError
        
        ## Interpolate function
        vector[V] interpolate(vector[T] x) except +ValueError
        
        ## Interpolate function
        V interpolate(T x) except +ValueError