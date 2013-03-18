/* wripaca.i */

%module wripaca 
%pythoncode %{
__author__="Jon Wright <wright@esrf.fr>"
__date__="March 2013"
__version__="?"
__doc__="..."


%}


%{
#include "../include/wripaca.h"
#include "../include/affine.h"
%}


%include "carrays.i"
%array_functions( double, doubleArray);


typedef double real;
typedef real rmatrix[9];
typedef real vector[3];
typedef struct {
    rmatrix m;
    vector  v;
} afmatrix;

typedef float scalar;


scalar determinant3( rmatrix );


/*
int omegacalc( vector , rmatrix , rmatrix , rmatrix ,
        vector , real , real* , real* );
*/
