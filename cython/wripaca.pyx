

cimport cwripaca

cimport numpy as np
import numpy as np


def determinant3( np.ndarray[np.double_t, ndim=1] arg ):
    return cwripaca.determinant3( <cwripaca.rmatrix> arg.data )



def omegacalc( np.ndarray[np.double_t, ndim=1] gve , 
        np.ndarray[np.double_t, ndim=1] pre , 
        np.ndarray[np.double_t, ndim=1] post ,
        np.ndarray[np.double_t, ndim=1] axis , 
        double wvln ):
    cdef cwripaca.real om1, om2
    success = cwripaca.omegacalc( 
            <cwripaca.vector> np.PyArray_DATA(gve), 
            <cwripaca.rmatrix> np.PyArray_DATA(pre), 
            <cwripaca.rmatrix> np.PyArray_DATA(post),
            <cwripaca.vector> np.PyArray_DATA(axis), 
            <cwripaca.real> wvln , 
            &om1, 
            &om2 )
    if success != 0:
        raise Exception("Badness in omegacalc %d"%(success))
    return om1, om2



