

cimport cwripaca

cimport numpy as np
import numpy as np


def determinant3( np.ndarray[np.double_t, ndim=1] arg ):
    return cwripaca.determinant3( <cwripaca.rmatrix> arg.data )



def omegacalc( np.ndarray[np.double_t, ndim=1, mode='c'] gve , 
        np.ndarray[np.double_t, ndim=1, mode='c'] pre , 
        np.ndarray[np.double_t, ndim=1, mode='c'] post ,
        np.ndarray[np.double_t, ndim=1, mode='c'] axis , 
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

def hklcalc_many( 
        np.ndarray[np.double_t, ndim=2, mode='c'] XL, # Diffraction spot in lab
        np.ndarray[np.double_t, ndim=1, mode='c'] axis, # rotation axis direction
        np.ndarray[np.double_t, ndim=1, mode='c'] angle, # omegas
        np.ndarray[np.double_t, ndim=1, mode='c'] UBI, # for h = UBI.g
        np.ndarray[np.double_t, ndim=1, mode='c'] T, # grain translation
        np.ndarray[np.double_t, ndim=2, mode='c'] hkl, # output
        double wvln
        ):
    cdef int success , i, j
    cdef cwripaca.vector x,dLn, O, k, g, h 
    cdef cwripaca.real M
    assert XL.shape[1] == 3
    for i in range( XL.shape[0] ):
        x[0] = XL[i,0]  
        x[1] = XL[i,1]
        x[2] = XL[i,2]
        cwripaca.hklcalc( x, 
                <cwripaca.vector*> np.PyArray_DATA(axis),
                <cwripaca.real> angle[i],
                <cwripaca.rmatrix> np.PyArray_DATA( UBI ),
                <cwripaca.vector> np.PyArray_DATA( T ),
                &dLn, &O, &M, &k, &g, &h, wvln)
        hkl[i,0] = h[0]
        hkl[i,1] = h[0]
        hkl[i,2] = h[0]


