
cimport cython
cimport cwripaca

cimport numpy as np
import numpy as np

cdef extern from "math.h":
    double floor(double x)
    double sqrt(double x)
    double atan2(double x, double y)
    double sin(double x)
    double cos(double x)

def determinant3( np.ndarray[np.double_t, ndim=1, mode='c'] arg ):
    return cwripaca.determinant3( <cwripaca.rmatrix> arg.data )


def rotate_vector_axis_angle(
        np.ndarray[np.double_t, ndim=1, mode='c'] axis,
        double angle,
        np.ndarray[np.double_t, ndim=1, mode='c'] inp,
        np.ndarray[np.double_t, ndim=1, mode='c'] outp,
        ):
   cwripaca.rotate_vector_axis_angle( 
           <cwripaca.vector> np.PyArray_DATA(axis),
           angle,
           <cwripaca.vector> np.PyArray_DATA(inp),
           <cwripaca.vector*> np.PyArray_DATA(outp),
        )   

def rotate_vectors_axis_angle(
        np.ndarray[np.double_t, ndim=1, mode='c'] axis,
        np.ndarray[np.double_t, ndim=1, mode='c'] angles,
        np.ndarray[np.double_t, ndim=2, mode='c'] inp,
        np.ndarray[np.double_t, ndim=2, mode='c'] outp,
        ):
    cdef int i
    assert angles.shape[0] == inp.shape[0]
    assert angles.shape[0] == outp.shape[0]
    assert inp.shape[1] == 3
    assert outp.shape[1] == 3

    for i in range(angles.shape[0]):
       cwripaca.rotate_vector_axis_angle( 
           <cwripaca.vector> np.PyArray_DATA(axis),
           <cwripaca.real>   angles[i],
           <cwripaca.vector> np.PyArray_DATA(inp[i]),
           <cwripaca.vector*> np.PyArray_DATA(outp[i]),
        )   



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

def anglediff(double a1, double a2):
    return atan2(sin(a1-a2), cos(a1-a2))


def omegacalcclose( 
        np.ndarray[np.double_t, ndim=2, mode='c'] gve , 
        np.ndarray[np.double_t, ndim=1, mode='c'] pre , 
        np.ndarray[np.double_t, ndim=1, mode='c'] post ,
        np.ndarray[np.double_t, ndim=1, mode='c'] axis , 
        np.ndarray[np.double_t, ndim=1, mode='c'] omegaobs,
        np.ndarray[np.double_t, ndim=1, mode='c'] omegacalc,
        np.ndarray[np.double_t, ndim=1, mode='c'] omegaerr,
        double wvln ):
    cdef cwripaca.real om1, om2, e1, e2
    cdef int i
    assert gve.shape[1] == 3
    assert gve.shape[0] == omegaobs.shape[0]
    assert gve.shape[0] == omegacalc.shape[0]
    assert gve.shape[0] == omegaerr.shape[0]
    for i in range(omegaobs.shape[0]):
        success = cwripaca.omegacalc( 
            <cwripaca.vector> np.PyArray_DATA(gve[i]), 
            <cwripaca.rmatrix> np.PyArray_DATA(pre), 
            <cwripaca.rmatrix> np.PyArray_DATA(post),
            <cwripaca.vector> np.PyArray_DATA(axis), 
            <cwripaca.real> wvln , 
            &om1, 
            &om2 )
        if success != 0:
            raise Exception("Badness in omegacalc %d"%(success))
        e1 = atan2( sin( om1 - omegaobs[i]), cos( om1 - omegaobs[i] ))
        e2 = atan2( sin( om2 - omegaobs[i]), cos( om2 - omegaobs[i] ))
        if abs(e1) < abs(e2):
            omegacalc[i] = om1
            omegaerr[i] = e1
        else:
            omegacalc[i] = om2
            omegaerr[i] = e2
        #        print i, omegaobs[i],om1,om2,e1,e2




def hklcalc_many( 
        np.ndarray[np.double_t, ndim=2, mode='c'] XL, # Diffraction spot in lab
        np.ndarray[np.double_t, ndim=1, mode='c'] axis, # rotation axis direction
        np.ndarray[np.double_t, ndim=1, mode='c'] angle, # omegas
        np.ndarray[np.double_t, ndim=1, mode='c'] UBI, # for h = UBI.g
        np.ndarray[np.double_t, ndim=1, mode='c'] T, # grain translation
        np.ndarray[np.double_t, ndim=1, mode='c'] pre, # 9x1 wedgechi
        np.ndarray[np.double_t, ndim=1, mode='c'] post, # 9x1 wedgechi
        np.ndarray[np.double_t, ndim=2, mode='c'] hkl, # output
        double wvln,
        np.ndarray[np.double_t, ndim=2, mode='c'] kcalc, # output
        ):
    cdef int success , i, j
    cdef cwripaca.vector x,dLn, O, k, g, h 
    cdef cwripaca.real M
    assert XL.shape[1] == hkl.shape[1] == kcalc.shape[1] == 3
    for i in range( XL.shape[0] ):
        x[0] = XL[i,0]  
        x[1] = XL[i,1]
        x[2] = XL[i,2]
        cwripaca.hklcalc( x, 
                <cwripaca.vector*> np.PyArray_DATA(axis),
                <cwripaca.real> angle[i],
                <cwripaca.rmatrix> np.PyArray_DATA( UBI ),
                <cwripaca.vector> np.PyArray_DATA( T ),
                <cwripaca.rmatrix> np.PyArray_DATA( pre ),
                <cwripaca.rmatrix> np.PyArray_DATA( post ),
                &dLn, &O, &M, 
                <cwripaca.vector*> np.PyArray_DATA( kcalc[i] ), 
                &g, &h, wvln)
        hkl[i,0] = h[0]
        hkl[i,1] = h[1]
        hkl[i,2] = h[2]



@cython.boundscheck(False)
@cython.wraparound(False)
def ih_drlv( np.ndarray[np.double_t, ndim=2, mode='c'] hkl ,
     np.ndarray[np.double_t, ndim=1, mode='c'] drlv ):
    """
    sets hkl to the nearest integer values 
    sets drlv to |hkl - hinteger|
    """
    cdef int i
    cdef double ih0, ih1, ih2
    assert drlv.shape[0] == hkl.shape[0]
    assert hkl.shape[1] == 3
    for i in range( hkl.shape[0] ):
        ih0 = floor(hkl[i,0]+0.5)
        ih1 = floor(hkl[i,1]+0.5)
        ih2 = floor(hkl[i,2]+0.5)
        drlv[i] = sqrt( (hkl[i,0]-ih0)*(hkl[i,0]-ih0) +
                        (hkl[i,1]-ih1)*(hkl[i,1]-ih1) +
                        (hkl[i,2]-ih2)*(hkl[i,2]-ih2) )
        hkl[i,0] = ih0
        hkl[i,1] = ih1
        hkl[i,2] = ih2
        
