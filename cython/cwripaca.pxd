


cdef extern from "../include/affine.h":
    ctypedef double real
    ctypedef real rmatrix[9]
    ctypedef real vector[3]
    ctypedef struct afmatrix:
        rmatrix m
        vector  v

    ctypedef real scalar

    scalar determinant3( rmatrix )

    int rotate_vector_axis_angle(vector, scalar, vector, vector*)

    int omegacalc( vector , rmatrix , rmatrix ,
        vector , real , real* , real* )

    int hklcalc( vector , vector *, real ,
      rmatrix , vector ,
      rmatrix, rmatrix,
      vector *, vector *, real *,
      vector *,
      vector *, 
      vector *,
      real )


