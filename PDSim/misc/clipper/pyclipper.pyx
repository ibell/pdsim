
import sys
import struct
from ctypes import *
import copy
import unicodedata
import time
import cython

from cython.operator cimport dereference as deref

cdef extern from "Python.h":
    Py_INCREF(object o)
    object Py_BuildValue(char *format, ...)
    object PyBuffer_FromMemory(void *ptr, int size)
    #int PyArg_ParseTuple(object struct,void* ptr)
    char* PyString_AsString(object string)
    int PyArg_VaParse(object args,  char *format, ...)
    int PyArg_Parse(object args,  char *format, ...)
    int PyObject_AsReadBuffer(object obj,  void* buffer, int* buffer_len)
    object PyBuffer_FromObject(object base, int offset, int size)
    object PyBuffer_FromReadWriteObject(object base, int offset, int size)
    PyBuffer_New(object o)

cdef extern from "stdio.h":
    cdef void printf(char*,...)

cdef extern from "stdlib.h":
    cdef void* malloc(unsigned int size)
    cdef void* free(void* p)
    char *strdup(char *str)
    int strcpy(void* str, void* src)
    int memcpy(void* str, void* src, int size)

from libcpp.vector cimport vector
        
cdef extern from "clipper.hpp" namespace "ClipperLib":
    #enum ClipType { ctIntersection, ctUnion, ctDifference, ctXor };

    cdef enum ClipType:
        ctIntersection= 1,
        ctUnion=2,
        ctDifference=3,
        ctXor=4
    
    #enum PolyType { ptSubject, ptClip };
    cdef enum PolyType:
        ptSubject= 1,
        ptClip=2
    
    # By far the most widely used winding rules for polygon filling are
    # EvenOdd & NonZero (GDI, GDI+, XLib, OpenGL, Cairo, AGG, Quartz, SVG, Gr32)
    # Others rules include Positive, Negative and ABS_GTR_EQ_TWO (only in OpenGL)
    # see http://glprogramming.com/red/chapter11.html
    #enum PolyFillType { pftEvenOdd, pftNonZero, pftPositive, pftNegative };
    cdef enum PolyFillType:
        pftEvenOdd= 1,
        pftNonZero=2,
        pftPositive=3,
        pftNegative=4
    
    ctypedef signed long long long64
    ctypedef unsigned long long ulong64
    ctypedef char bool 

    cdef struct IntPoint:
        long64 X
        long64 Y
        #IntPoint(X, Y)
        #IntPoint(long64 X = 0, long64 Y = 0)
        #IntPoint(long64 x = 0, long64 y = 0): X(x), Y(y) {};
        #friend std::ostream& operator <<(std::ostream &s, IntPoint &p);

    cdef cppclass Polygon:
        Polygon()
        void push_back(IntPoint&)
        IntPoint& operator[](int)
        IntPoint& at(int)
        int size()
        
    cdef cppclass Polygons:
        Polygons()
        void push_back(Polygon&)
        Polygon& operator[](int)
        Polygon& at(int)
        int size()

    #ctypedef std::vector< IntPoint > Polygon
    #ctypedef std::vector< Polygon > Polygons
    #ctypedef vector[IntPoint] Polygon
    #ctypedef vector[Polygon] Polygons
    #ctypedef  IntPoint Polygon
    #ctypedef  Polygon Polygons

    #std::ostream& operator <<(std::ostream &s, Polygon &p);
    #std::ostream& operator <<(std::ostream &s, Polygons &p);

    cdef struct ExPolygon:
        vector[IntPoint] outer
        vector[IntPoint] holes

    #ctypedef std::vector< ExPolygon > ExPolygons
    cdef cppclass ExPolygons:
        ExPolygons()
        void push_back(Polygon&)
        ExPolygon& operator[](int)
        ExPolygon& at(int)
        int size()  
          
    #ctypedef vector[ExPolygon] ExPolygons
    #ctypedef ExPolygon ExPolygons

    #enum JoinType { jtSquare, jtRound, jtMiter };
    cdef enum JoinType:
        jtSquare= 1,
        jtRound=2
        jtMiter=3
        
    bool ReverseSolution(bool value)
    
    #bool Orientation(const Polygon &poly);
    bool Orientation(vector[IntPoint] poly)

    #double Area(const Polygon &poly)
    double Area(Polygon poly)

    #void OffsetPolygons(const Polygons &in_polys, Polygons &out_polys, double delta, JoinType jointype = jtSquare, double MiterLimit = 2);
    void OffsetPolygons(Polygons in_polys, Polygons out_polys, double delta, JoinType jointype, double MiterLimit)
    
    #void SimplifyPolygon(const Polygon &in_poly, Polygons &out_polys);
    void SimplifyPolygon(Polygon in_poly, Polygons out_polys)

    #void SimplifyPolygons(const Polygons &in_polys, Polygons &out_polys);
    void SimplifyPolygons(Polygons in_polys, Polygons out_polys)

    #void SimplifyPolygons(Polygons &polys);
    void SimplifyPolygons(Polygons polys)

    #void ReversePoints(Polygon& p);
    void ReversePoints(Polygon p)

    #void ReversePoints(Polygons& p);
    void ReversePoints(Polygons p)

    #used internally ...

    #enum EdgeSide { esNeither = 0, esLeft = 1, esRight = 2, esBoth = 3 };
    cdef enum EdgeSide:
        esNeither=0,
        esLeft= 1,
        esRight=2,
        esBoth=3

    #enum IntersectProtects { ipNone = 0, ipLeft = 1, ipRight = 2, ipBoth = 3 };
    cdef enum IntersectProtects:
        ipNone=0,
        ipLeft= 1,
        ipRight=2,
        ipBoth=3

    struct TEdge:
        long64 xbot
        long64 ybot
        long64 xcurr
        long64 ycurr
        long64 xtop
        long64 ytop
        double dx
        long64 tmpX
        PolyType polyType
        EdgeSide side 
        int windDelta #1 or -1 depending on winding direction
        int windCnt
        int windCnt2 #winding count of the opposite polytype
        int outIdx
        TEdge *next
        TEdge *prev
        TEdge *nextInLML
        TEdge *nextInAEL
        TEdge *prevInAEL
        TEdge *nextInSEL
        TEdge *prevInSEL

    struct IntersectNode:
        TEdge          *edge1
        TEdge          *edge2
        IntPoint        pt
        IntersectNode  *next

    struct LocalMinima:
        long64        Y
        TEdge        *leftBound
        TEdge        *rightBound
        LocalMinima  *next

    cdef struct Scanbeam:
        long64    Y
        Scanbeam *next

    cdef struct OutRec:
        int     idx
        bool    isHole
        OutRec *FirstLeft
        OutRec *AppendLink
        OutPt  *pts
        OutPt  *bottomPt
        OutPt  *bottomFlag
        EdgeSide sides

    cdef struct OutPt:
        int     idx
        IntPoint pt
        OutPt   *next
        OutPt   *prev

    cdef struct JoinRec:
        IntPoint  pt1a
        IntPoint  pt1b
        int       poly1Idx
        IntPoint  pt2a
        IntPoint  pt2b
        int       poly2Idx

    cdef struct HorzJoinRec:
        TEdge    *edge
        int       savedIdx



    #ctypedef std::vector < OutRec* > PolyOutList
    #ctypedef  OutRec PolyOutList

    #ctypedef std::vector < TEdge* > EdgeList
    #ctypedef TEdge EdgeList 

    #ctypedef std::vector < JoinRec* > JoinList
    #ctypedef std::vector < HorzJoinRec* > HorzJoinList
    #ctypedef JoinRec JoinList
    #ctypedef HorzJoinRec HorzJoinList

    cdef struct IntRect:
        long64 left
        long64 top
        long64 right
        long64 bottom
        
    cdef cppclass Clipper:
        Clipper()
        #~Clipper()
        bool Run(ClipType clipType, Polygons solution, PolyFillType subjFillType, PolyFillType clipFillType)
        bool Execute(ClipType clipType,  ExPolygons solution,  PolyFillType subjFillType,  PolyFillType clipFillType)
        void Clear()
        bool ReverseSolution()
        void ReverseSolution(bool value)

        # Inherited methods from Clipper Base
        #bool AddPolygon(const Polygon &pg, PolyType polyType)
        #bool AddPolygons( const Polygons &ppg, PolyType polyType)

        bool AddPolygon( Polygon pg, PolyType polyType)
        bool AddPolygons( Polygons ppg, PolyType polyType)
        void Clear()
        IntRect GetBounds()

#=============================  Namespace methods ==================
INTERSECTION = ctIntersection
DIFFERENCE = ctDifference
UNION = ctUnion
XOR = ctXor

#===========================================================

#===========================================================
# OffsetPolygons(const Polygons &in_polys, Polygons &out_polys,  double delta, JoinType jointype = jtSquare, double MiterLimit = 2);
def offset( pypolygons, delta=100,  jointype = jtSquare, double MiterLimit = 2):
    print("Offset polygon")
    cdef Polygon poly =  Polygon() 
    cdef IntPoint a
    cdef Polygons polys =  Polygons()
    for pypolygon in pypolygons:
        for pypoint in pypolygon:
            a = IntPoint(pypoint[0], pypoint[1])
            poly.push_back(a)
     
        polys.push_back(poly)  

    cdef Polygons solution
    OffsetPolygons( polys, solution,  delta,  jointype, MiterLimit)
    n = solution.size()
    sol = []
    print("Solution is made of %i loops"%n )

    cdef IntPoint point
    for i in range(n):
        poly = solution[i]
        m = poly.size()
        print("loop has %i points"%m)
        loop = []
        for i in range(m):
            point = poly[i]
            print(point.X ,point.Y)
            loop.append([point.X ,point.Y])
        sol.append(loop)
    return sol

#===========================================================
# void SimplifyPolygons(const Polygons &in_polys, Polygons &out_polys)
def simplify_polygons(pypolygons):
    cdef Polygon poly =  Polygon() 
    cdef IntPoint a
    cdef Polygons polys =  Polygons()
    for pypolygon in pypolygons:
        for pypoint in pypolygon:
            a = IntPoint(pypoint[0], pypoint[1])
            poly.push_back(a)
        polys.push_back(poly)  

    cdef Polygons solution
    SimplifyPolygons( polys, solution)
    n = solution.size()
    sol = []

    cdef IntPoint point
    for i in range(n):
        poly = solution[i]
        m = poly.size()
        loop = []
        for i in range(m):
            point = poly[i]
            loop.append([point.X ,point.Y])
        sol.append(loop)
    return sol

#=============================  Pyclipper ==================
cpdef double area(pypolygon):

    cdef Polygon square =  Polygon() 
    cdef IntPoint a
    for p in pypolygon:
        a = IntPoint(p[0], p[1])
        square.push_back(a)
    return Area(square)
        
#===========================================================
cdef class Pyclipper:
    cdef Clipper *thisptr      # hold a C++ instance which we're wrapping
    error_code = {-1:"UNSPECIFIED_ERROR", -2: "BAD_TRI_INDEX", -3:"NO_VOX_MAP", -4:"QUERY_FAILED"}

    #===========================================================
    def __cinit__(self):
        self.thisptr = new Clipper()

    #===========================================================
    def __dealloc__(self):
        del self.thisptr
        

    #===========================================================
    #bool AddPolygon(Polygon pg, PolyType polyType)
    def subject_polygon(self, pypolygon):

        cdef Polygon square =  Polygon() 
        cdef IntPoint a
        for p in pypolygon:
            a = IntPoint(p[0], p[1])
            square.push_back(a)

        cdef Polygons subj =  Polygons() 
        subj.push_back(square)  
        self.thisptr.AddPolygons(subj, ptSubject)


    #===========================================================
    #bool AddPolygon(Polygon pg, PolyType polyType)
    def clip_polygon(self, pypolygon):
        cdef Polygon square =  Polygon() 
        cdef IntPoint a
        for p in pypolygon:
            a = IntPoint(p[0], p[1])
            square.push_back(a)

        cdef Polygons subj =  Polygons() 
        subj.push_back(square)  
        self.thisptr.AddPolygons(subj, ptClip)
        
    #===========================================================
    # IntRect GetBounds();
    def GetBounds(self):
        cdef IntRect rect=  self.thisptr.GetBounds()
        return {"top":rect.top, "left":rect.left, "right":rect.right, "bottom":rect.bottom}
    
    #===========================================================
    # IntRect GetBounds();
    def Clear(self):
        #self.thisptr.Clear()
        return 
        
    #===========================================================
    # reverse solution
    def ReverseSolution(self, direction):
        #cdef bool ret = self.thisptr.ReverseSolution()
        return

    #===========================================================
    # Execute(ClipType clipType, vector[vector[IntPoint]] solution, PolyFillType subjFillType, PolyFillType clipFillType)
    def execute(self, mode=ctDifference):
        cdef Polygons solution
        self.thisptr.Run(mode, solution, pftEvenOdd , pftEvenOdd)
        n = solution.size()
        sol = []
        cdef Polygon poly
        cdef IntPoint point
        for i in range(n):
            poly = solution[i]
            m = poly.size()
            loop = []
            for i in range(m):
                point = poly[i]
                loop.append([point.X ,point.Y])
            sol.append(loop)
        return sol
        
    #===========================================================
    # test()
    def test(self):
        # Import polygons to structures
        #read_polygons_MATLAB(prhs[0], subj);
        #read_polygons_MATLAB(prhs[1], clip);
        
        cdef Polygon hole =  Polygon() 
        cdef IntPoint b
        for p in [[-1,-1], [6, 5], [6 , 6], [5, 6]]:
            b = IntPoint(p[0], p[1])
            hole.push_back(b)
            
        print(hole.size())

        cdef Polygons clip =  Polygons() 
        clip.push_back(hole)  
        self.thisptr.AddPolygons(clip, ptClip)
        
        
        cdef Polygons solution

        self.thisptr.Run(ctDifference, solution, pftEvenOdd, pftEvenOdd)
        n = solution.size()
        print("Solution is made of %i loops"%n  )
        cdef Polygon poly
        cdef IntPoint point
        for i in range(n):
            poly = solution[i]
            m = poly.size()
            print("First loop has %i points"%m)
            
            for i in range(m):
                point = poly[i]
                print(point.X ,point.Y)



