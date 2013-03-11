#include "affine.h"
#include <stdio.h>

void printvec3( char* name, real vec[] ){
    printf("%s %f %f %f\n",name,vec[0],vec[1],vec[2]);
}

/**
 * Given XL, YL, ZL, Omega, UBI, T, find h,k,l
 *
 * hkl = UBI.gvector
 */
int hklcalc( vector XL, vector *axis, real ang,
            rmatrix UBI, vector T,
      vector *dLn, vector *O, real *M,
      vector *kvector,
      vector *gvector, 
      vector *hkl,
      real wvln ){
    int err;

    err = rotate_vector_axis_angle( axis, -ang, T, O ); 
    if ( err != 0 ) return err;

#ifdef DEBUG
    printvec3("axis",*axis);
    printvec3("T",T);
    printvec3("O",*O);
    printvec3("XL",XL);
#endif
    /* Unit vector in lab is: */
    vec3sub( XL, *O, *dLn );
#ifdef DEBUG
    printvec3("dLn",*dLn);
#endif
    *M = _norm3(*dLn);
#ifdef DEBUG
    printf("M %f   wvln %f\n",*M,wvln);
#endif
    vec3sdiv( *dLn, *M, *dLn );
#ifdef DEBUG
    printvec3("norm dLn",*dLn);
#endif
    *dLn[0] = *dLn[0] - 1.0; /* Beam along x */
    vec3sdiv( *dLn, wvln, *kvector );    
#ifdef DEBUG
    printvec3("kvec",*kvector);
#endif

    err = rotate_vector_axis_angle( axis, ang, *kvector, gvector );
    if ( err != 0 ) return err;
    
    matVec( UBI, *gvector, *hkl );
#ifdef DEBUG
    printvec3("gvec",*gvector);
#endif

    return 0;
}

/**
 * Determine the omega angle for reflection h,k,l given UBI, pre, post and
 * current rotation axis
 *
 * From ImageD11/gv_general.py
 */
int omegacalc( vector hkl, rmatrix UBI, rmatrix pre, rmatrix post,
        vector axis, real wvln, real *ang1, real *ang2){
/**
 *     Get the k and angles given the g in 
 *        g = pre . rot(axis, angle) . post . k
 *   ...where k will satisfy the Laue equations
 *   
 *   The recipe is:
 *       Find the components of g along and perpendicular to the rotation axis
 *           co-ords are a0 = axis, 
 *                       a1 = axis x g, 
 *                       a2 = a1 x axis = ( axis x g ) x axis
 *       Rotated vector will be a0 + a1 sin(angle) + a2 cos(angle)
 *       Laue condition says that [incident beam].[k] = k.sin(theta)
 *                                                    = 2 sin^2(theta) / lambda 
 *                                                    = sin(t)/d = |g|.sin(t)
 *                                                    = |g|*|g|*lambda / 2
 *       Apply any post rotations to the incident beam 
 *       Apply any pre-rotations to the g vector
 *       |g| = [a0 + a1 sin(angle) + a2 cos(angle) ] . [incident beam]
 *      => solve for angle
 *
 *       http://www.ping.be/~ping1339/gonio.htm
 *       a sin(u) + b cos(u) = c
 *       let: tan(u') = -b/a
 *       and: A = a / cos(u')
 *       Finally you'll find:
 *            A.sin(u-u') = c
 *            A.cos(pi/2 - u - u') = c
 */
    vector tmp, rotg, rbeam;
    vector a0, a1, a2;
    real rbda0, rbda1, rbda2, modg, kdotbeam, p, k, x_plus_p;
    real q; /* Expect these to be optimised away... */
    /* Step 1: compute the g-vector */
    matVec( UBI, hkl, rotg );
    /* Having this as inline gives a 100X speedup - it must be 
     * optimising something away */
    modg = norm3( rotg );
    kdotbeam = -modg*modg/2.;
    /* Apply the pre-rotation */
//    mat3_transform_vec( pre, tmp, &rotg );
    /* Put incident beam in tmp */
    rbeam[0] = -1.0/wvln;
    rbeam[1] = 0.0;
    rbeam[2] = 0.0;
    /* Rotate incident beam (inverse is transpose for rotation) */
//    mat3_T_transform_vec( post, tmp, &rbeam);
    /* Now find the components of g w.r.t our rotation axis */
    /* a1 = perpendicular to both axis and g */
    crossProduct( axis, rotg , a1 );
    /* a2 = perpendicular to axis and along g */
    crossProduct( a1, axis, a2 );
    vec3sub( rotg, a2, a0 ); 
    /* Dot products of these 3 with the incident beam */
    rbda0 = _dot3( rbeam, a0 );
    rbda1 = _dot3( rbeam, a1 );
    rbda2 = _dot3( rbeam, a2 );
    /* k.b = rbda0 + rbda1.sin(t) + rbda2.cos(t) */
    /* a = rbda1;
       b = rbda2;
       c = kdotbeam - rbda0; */
    /* From wikipedia: 
    # http://en.wikipedia.org/wiki/List_of_trigonometric_identities#Linear_combinations
    # a.sin(x) + b.cos(x) = sqrt(a*a+b*b) sin(x+p)
    # with p = atan(b/a) */
    p = atan2( rbda2, rbda1 );
    k = sqrt( rbda1*rbda1 + rbda2*rbda2 );
    q = (kdotbeam - rbda0)/k;
    if ( q > 1 ) return 1;  /* No solution exists */
    if ( q < -1 ) return 1;
    x_plus_p = asin( q );
    *ang1 = x_plus_p + p;
    *ang2 = M_PI - x_plus_p + p;
    return 0;
}


#include "windows.h"
/*BOOL WINAPI QueryPerformanceCounter(
  _Out_  LARGE_INTEGER *lpPerformanceCount
);*/


int main(){
    /* inputs */
    vector XL={100.,100.,100.};
    vector axis = {0.,0.,1.};
    real ang = 0.3;
    rmatrix UBI = { 1.,0.,0.,0.,1.,0.,0.,0.,1.};
    vector T = {10,11,12};
    real wvln = 0.154;
    /* outputs */
    vector kvector, gvector, hkl, dL, O;
    rmatrix UB;
    real  M, ang1, ang2;
    vector kcalc, gcalc,dLcalc,Ocalc,Tcalc;
    vector ih;
    int err, ihkl[3];
    
    LARGE_INTEGER start, end, freq;

    printf("Hello\n");

    /* 7 seconds for 2000 grains and 60000 peaks */
    QueryPerformanceCounter(&start);
    for(err=0;err<60000*2000;err++)
      hklcalc(XL, &axis, ang, UBI, T, &dL, &O, &M, &kvector, &gvector, &hkl, wvln );
    QueryPerformanceCounter(&end);
    QueryPerformanceFrequency(&freq);
    printf("Cycles %d %d %f /s\n",end.LowPart-start.LowPart, freq.LowPart,
            (double)(end.LowPart-start.LowPart)/(double)freq.LowPart);
    printvec3("dL",dL);
    printf("M %f\n",M);
    printf("hkl %f %f %f ",hkl[0],hkl[1],hkl[2]);
    round3( hkl, ih );
    printf("hkl %f %f %f \n",ih[0],ih[1],ih[2]);
    
    err = inverse_mat3( UBI, &UB );
    if (err != 0) return err;

/*    scanf("%d",&(ihkl[0]));
    scanf("%d",&(ihkl[1]));
    scanf("%d",&(ihkl[2]));
    ih[0] = (double) ihkl[0];
    ih[1] = (double) ihkl[1];
    ih[2] = (double) ihkl[2];*/
    /* Compute the ideal angles given hkl here */
    QueryPerformanceCounter(&start);
    /* 10.5 seconds for 2000 grains and 60000 peaks */
    for(err=0;err<60000*2000;err++){
     ih[0]=(real) (err%10);
     ih[1]=(real) (err%11);
     ih[2]=(real) (err%12);
     omegacalc(ih, UB, UB, UB,
        axis, wvln,  &ang1, &ang2);
    }
    QueryPerformanceCounter(&end);
    QueryPerformanceFrequency(&freq);
    printf("Computed for (%f, %f, %f)\n",ih[0],ih[1],ih[2]);
    printf("Cycles %d %d %f /s\n",end.LowPart-start.LowPart, freq.LowPart,
            (double)(end.LowPart-start.LowPart)/(double)freq.LowPart);
    printf("Ideal angles are %f %f radians\n",ang1, ang2);
    printf("Ideal angles are %f %f degrees\n",ang1*180.0/M_PI, ang2*180.0/M_PI);


    return 0;
}
