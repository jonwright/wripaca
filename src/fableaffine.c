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
    for(err=0;err<60000*200;err++)
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
    for(err=0;err<60000*200;err++){
     ih[0]=(real) (err%10);
     ih[1]=(real) (err%11);
     ih[2]=(real) (err%12);
     omegacalc(ih, UB, UB,
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
