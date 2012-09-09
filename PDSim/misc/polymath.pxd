import cython
cimport cython

cimport numpy as np
import numpy as np

from cpython cimport array
import array

cdef class LineSegment:
    cdef double x0,x1,y0,y1,y2,y3,x_i,y_i
    
    cpdef bint intersects(self, LineSegment LS2)
    
cdef class Polygon:
    cdef public np.ndarray x, y
    cdef public list xy
    
    cpdef bint inpoly(self, double x, double y)
    cpdef double area(self)

cdef class PolygonOperator: 
    cdef Polygon poly1, poly2
    cdef list Nodes, Parts, in1, in2
    cdef public dict new_poly, unmatched_partners
    
    cpdef list AND(self)
    cpdef list OR(self)
    cpdef list XOR(self)
    cpdef intersect(self)
    cpdef make_parts(self, Polygon poly, list sorted)
    cpdef attach_nodes(self)
    cpdef get_partners(self, PolyPart part)
    cpdef more_partners(self)
    cpdef match_parts(self)
    cpdef remove_dup_poly(self, dict polydict)
    
cdef class Node:
    cdef public double x, y
    cdef public long I
    cdef public list children
    
cdef class PolyPart:
    cdef public np.ndarray x,y
    cdef public Polygon parent
    cdef public tuple nodes
    
    cpdef bint is_within(self, Polygon poly)