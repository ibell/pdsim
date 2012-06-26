%module GeometryCore
%{
#include "GeometryCore.h"
%}

%include "typemaps.i"

//All of the parameters that match the list below will be converted to output parameters
%apply double *OUTPUT { double *nx, double *ny, double *x, double *y, double *V, double *dV, double *cx, double *cy, double *A, int *errCode };

%include "GeometryCore.h"