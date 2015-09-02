from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cython
import itertools

if (matplotlib.__version__ < '1.2'):
    from matplotlib.nxutils import pnpoly
else:
    from matplotlib.path import Path as mpl_path




"""
A module for doing polygon math - finding regions of overlap and differences of
two polygons

Updated Davide Ziviani 9/2015 
- The matplotlib.nxutils module has been removed from Matplotlib. Use "matplotlib.path.Path.contains_point" instead
  you make a Path object using the vertices, and then use its contains_point() method:
  
  path = Path(polygonVerts)
  isInside = path.contains_point(point) 


Info on pnpoly:
pnpoly(...)
inside = pnpoly(x, y, xyverts)
 
Return 1 if x,y is inside the polygon, 0 otherwise.
 
*xyverts*
    a sequence of x,y vertices.

"""

cdef class Polygon(object):
    def __init__(self, np.ndarray x, np.ndarray y):
        self.x = x
        self.y = y
        self.xy = zip(self.x, self.y)
            
    cpdef bint inpoly(self, double testx, double testy):
        """
        
        Recommendation from 
        http://www.heikkitoivonen.net/blog/2009/01/26/point-in-polygon-in-python/
        """
        #TODO: to fix
        if(matplotlib.__version__ < '1.2'):
            return pnpoly(testx, testy, self.xy)
        
        return mpl_path(self.xy).contains_points([(testx,testy)])


    cpdef double area(self):
        cdef long i
        cdef double area = 0.0
        for i in range(len(self.x)-1):
            area += self.x[i]*self.y[i+1] - self.y[i]*self.x[i+1]
        return area/2.0

    def plot(self, ax, **kwargs):
        ax.plot(self.x, self.y, **kwargs)
    
    def fill(self, ax, **kwargs):
        ax.fill(self.x, self.y, alpha = 0.5, **kwargs)
    
cdef class LineSegment(object):
    """ A line segment with two nodes """
    @cython.locals(x0 = cython.double, 
                   y0 = cython.double,
                   x1 = cython.double, 
                   y1 = cython.double)
    def __init__(self,x0,y0,x1,y1):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        
        #Intersection point
        self.x_i = 1e99
        self.y_i = 1e99 
        
    @cython.cdivision(True)
    cpdef bint intersects(self, LineSegment LS2):
        cython.declare(x0=cython.double,
                       y0=cython.double,
                       x1=cython.double,
                       y1=cython.double,
                       x2=cython.double,
                       y2=cython.double,
                       x3=cython.double,
                       y3=cython.double
                       )
        x0 = self.x0
        x1 = self.x1
        y0 = self.y0
        y1 = self.y1
        x2 = LS2.x0
        x3 = LS2.x1
        y2 = LS2.y0
        y3 = LS2.y1
        
        #Using Cramer's rule
        # http://en.wikipedia.org/wiki/Cramer%27s_rule
        # Upper left element of A
        a = x1-x0
        # Upper right element of A
        b = -(x3-x2)
        # Lower left element of A
        c = y1-y0
        # Lower right element of A
        d = -(y3-y2)
        # Top RHS element
        e = -(x0-x2)
        # Bottom RHS element
        f = -(y0-y2)
        
        if (a*d-b*c == 0.0):
            return False
        
        #t parameter for the first line
        t = (e*d-b*f)/(a*d-b*c)
        #s parameter for the second line
        s = (a*f-e*c)/(a*d-b*c)
        
        if 0 <= t <= 1 and 0 <= s <= 1:
            self.x_i = (x1-x0)*t + x0
            self.y_i = (y1-y0)*t + y0
            return True
        else:
            return False
        
cdef class PolygonOperator(object):
    
    def __init__(self,Polygon poly1, Polygon poly2):
        self.poly1 = poly1
        self.poly2 = poly2
        self.Nodes = []
        self.Parts = []
        self.new_poly = {}
        self.unmatched_partners = {}
        self.in1 = []
        self.in2 = []
   
    cpdef list AND(self):
        """
        Return all the new polygons that are within both polygon1 and polygon2 
        """
        polys = self.intersect()
        return [poly for poly in polys if poly in self.in1 and poly in self.in2] 
        
        
    cpdef list OR(self):
        """
        Return all the new polygons that are within either polygon1 and polygon2 
        """
        polys = self.intersect()
        return [poly for poly in polys if poly in self.in1 or poly in self.in2]
    
    cpdef list XOR(self):
        """
        Return all the new polygons that are within one of polygon1 and polygon2 but not within the other one 
        """
        polys = self.intersect()
        return [poly for poly in polys if poly in self.in1 and not poly in self.in2 or poly not in self.in1 and poly in self.in2]
    
    cpdef list Only1(self):
        """
        Return all the new polygons that are only within polygon 1 and not within polygon 2 
        """
        polys = self.intersect()
        return [poly for poly in polys if poly in self.in1 and not poly in self.in2]
    
    cpdef list Only2(self):
        """
        Return all the new polygons that are only within polygon 2 and not within polygon 1 
        """
        polys = self.intersect()
        return [poly for poly in polys if poly in self.in2 and not poly in self.in1]
    
    cpdef intersect(self):
        """
        
        """
        cdef list intersections = []
        cdef long i2
        cdef double xp, yp
        cdef Polygon poly
        
        cython.declare(i = cython.long,
                       j = cython.long 
                       ) 
        
        # First find the intersections between the curves
        for i in range(len(self.poly1.x)-1):
            LS1 = LineSegment(self.poly1.x[i], self.poly1.y[i], self.poly1.x[i+1], self.poly1.y[i+1]) 
            for j in range(len(self.poly2.x)-1):
                LS2 = LineSegment(self.poly2.x[j], self.poly2.y[j], self.poly2.x[j+1], self.poly2.y[j+1])
                if LS1.intersects(LS2):
                    intersections.append((LS1.x_i,LS1.y_i,i,j))
        
        #No intersections
        if not intersections:
            return []
        
        #Add node entries for each node that splits a line segment
        for I,(x_i,y_i,i,j) in enumerate(intersections):
            self.Nodes.append(Node(x_i,y_i,I))
            
        #Replace the coordinates with the node
        intersections = [(i,j,_Node) for (x_i,y_i,i,j),_Node in zip(intersections,self.Nodes)]
        
        import operator
        # Sort the insertions by the index for poly1
        sorted1 = sorted(intersections, key=operator.itemgetter(0))
        # Sort the insertions by the index for poly2
        sorted2 = sorted(intersections, key=operator.itemgetter(1))
        
        # For the first polygon, split it into parts (lines)
        self.make_parts(self.poly1, sorted1)
        # For the second polygon, split it into parts (lines)
        self.make_parts(self.poly2, sorted2)
        
        #Attach each line segment to a node
        self.attach_nodes()

        #Match parts to form new polygons
        self.match_parts()
        
        #Figure out which new polygons are in which old polygon
        for poly in self.new_poly.values():
            #Get a point that is hopefully in the new polygon
            i2 = 0
            xp = (poly.x[0]+poly.x[i2])/2.0
            yp = (poly.y[0]+poly.y[i2])/2.0
            while not poly.inpoly(xp, yp):
                i2 += 1
                xp = (poly.x[0]+poly.x[i2])/2.0
                yp = (poly.y[0]+poly.y[i2])/2.0
                if i2 == len(poly.x)-1:
                    raise ValueError('could not obtain a point in the polygon')
            
            if self.poly1.inpoly(xp, yp):
                self.in1.append(poly)
            if self.poly2.inpoly(xp, yp):
                self.in2.append(poly)
        
        return self.new_poly.values()
        
    cpdef remove_dup_poly(self, dict polydict):
        #Get all the keys
        keys = polydict.keys()
        #Sort each key
        skeys = [sorted(_keys) for _keys in keys]
        
        dups = dict()
        for key in polydict:
            if skeys.count(sorted(key)) > 1:
                #It is a duplicate, save it in a dictionary
                #Dictionary used because there can be no repetition
                dups[tuple(sorted(key))] = None
        
        #Remove the duplicate polygons
        for dup in dups:
            polydict.pop(dup)
        
    cpdef more_partners(self):
        """
        Do the parts that have more than just two curves forming their outlines
        """
        # unmatched_partners is a dict with key of unmatched lines, 
        # and value a tuple of the neighboring lines
        # Start with one part, doesn't matter which
        if not self.unmatched_partners:
            return
        
        start_part = self.unmatched_partners.keys()[0]
        parts = [start_part]
        
        while len(parts) == 1 or not start_part == parts[-1]:
            #Get the next 2 elements
            partners = self.unmatched_partners[parts[-1]]
            #Pick the one that is not the previous element
            if len(parts) >= 2 and not partners[0] == parts[-2]:
                parts.append(partners[0])
            else:
                parts.append(partners[1])
            #Remove the part that was before this part
            self.unmatched_partners.pop(parts[-2])
        
        # Remove the end element
        parts.pop(len(parts)-1)
        
        #Join the curves together
        x = parts[0].x
        y = parts[0].y
        
        for i, part in enumerate(parts):
            if i>0:
                oldpart = parts[-1]
                if abs(x[-1]-part.x[0])>1e-12 or abs(y[-1]-part.y[0])>1e-12:
                    #Flip the partner matrices 
                    x = np.r_[x, part.x[::-1]]
                    y = np.r_[y, part.y[::-1]]
                else:
                    #Keep the same orientation
                    x = np.r_[x, part.x]
                    y = np.r_[y, part.y]
                
        self.new_poly[tuple(parts)] = Polygon(x,y)
    
    cpdef match_parts(self):
        
        for part in self.Parts:
            #Get its partners
            self.get_partners(part)
            
        #Match any of the other parts that remain
        self.more_partners()
           
        #Remove any polygons with the same lines traversed in different ways   
        self.remove_dup_poly(self.new_poly)
            
    @cython.nonecheck(False)
    cpdef get_partners(self, PolyPart part):
        """
        Get the partners on the other polygon
        """
        #: list of partners to remove
        cdef list del_list = []
        cdef list partners = []
        cdef PolyPart partner
        cdef Node node
        
        for node in part.nodes:
            for part2 in node.children:
                #They are on opposite polygons
                if not part2.parent == part.parent:
                    #See if other part is within the polygon this part comes from
                    if part2.is_within(part.parent):
                        partners.append(part2)
        
        for partner in partners:
            # If if has a partner that shows up twice, it is a closed polygon, add this polygon
            if partners.count(partner) == 2:
                #You are going to remove this when finished
                del_list.append(partner)
                if abs(part.x[-1]-partner.x[0])>1e-12 or abs(part.y[-1]-partner.y[0])>1e-12:
                    #Flip the partner matrices 
                    x = np.r_[part.x, partner.x[::-1]]
                    y = np.r_[part.y, partner.y[::-1]]
                else:
                    x = np.r_[part.x, partner.x]
                    y = np.r_[part.y, partner.y]    
                self.new_poly[(part,partner)] = Polygon(x,y)
        
        for partner in del_list:
            partners.remove(partner)
        
        if partners:
            self.unmatched_partners[part] = tuple(partners)
                
    cpdef attach_nodes(self):
        """ 
        Find the parts attached to each node 
        """
        
        for node in self.Nodes:
            for part in self.Parts:
                if node in part.nodes:
                    node.children.append(part)
            
    cpdef make_parts(self, Polygon poly, list sorted):
        cdef long I, parent
        
        if poly == self.poly1:
            parent = 1
        elif poly == self.poly2:
            parent = 2

        for I in range(0,len(sorted)-1):
            Node1 = sorted[I][2]
            Node2 = sorted[I+1][2]
            
            # The first index in the polygon part
            #
            # The value stored in sorted is the index to the left of the 
            # intersection, so +1 to get the value in the segment
            IL = sorted[I][parent - 1]+1
            # The last index in the polygon part is the one that is to the left
            # of the next intersection
            IR = sorted[I+1][parent-1]
            #Insert the intersection points into the arrays
            xx = np.r_[Node1.x, poly.x[IL:IR+1], Node2.x]
            yy = np.r_[Node1.y, poly.y[IL:IR+1], Node2.y]
            
            self.Parts.append( PolyPart(xx, yy, nodes = (Node1,Node2), parent = poly) )
            
        #Last part is the one that goes from the last node to end, and wraps around to first node again
         
        Node1 = sorted[-1][2]
        Node2 = sorted[0][2]
        
        # The first index in the polygon part
        #
        # The value stored in sorted is the index to the left of the 
        # intersection, so +1 to get the value in the segment
        IL = sorted[-1][parent - 1]+1
        # The last index in the polygon part is the one that is to the left
        # of the next intersection
        IR = sorted[0][parent-1]
        #Insert the intersection points into the arrays
        xx = np.r_[Node1.x, poly.x[IL::], poly.x[0:IR+1], Node2.x]
        yy = np.r_[Node1.y, poly.y[IL::], poly.y[0:IR+1], Node2.y]
        
        self.Parts.append( PolyPart(xx, yy, nodes = (Node1,Node2), parent = poly) )
        
cdef class Node(object):
    """
    An intersection point between two lines
    """
    
    def __init__(self, x, y, I):
        self.x = x
        self.y = y
        self.I = I
        self.children = []
    
    def plot(self, ax): 
        ax.plot(self.x,self.y,'o')
            
cdef class PolyPart(object):
    
    def __init__(self, np.ndarray x, np.ndarray y, tuple nodes, Polygon parent):
        self.x = x
        self.y = y
        self.nodes = nodes
        self.parent = parent
    
    cpdef bint is_within(self, Polygon poly):
        cdef long IM = len(self.y)//2
        
        if IM > 0:
            return poly.inpoly(self.x[IM], self.y[IM])
        else:
            return 0
    
    def plot(self, ax):
        ax.plot(self.x, self.y,'-x')
            
            
            
            
            
            
        