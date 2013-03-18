


cdef extern from "../include/affine.h":
    ctypedef double real
    ctypedef real rmatrix[9]
    ctypedef real vector[3]
    ctypedef struct afmatrix:
        rmatrix m
        vector  v

    ctypedef real scalar

    scalar determinant3( rmatrix )

    int omegacalc( vector , rmatrix , rmatrix ,
        vector , real , real* , real* )



