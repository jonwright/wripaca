#include "affine.h"
#include <stdio.h>

/**
 * Given XL, YL, ZL, Omega, UBI, T, find h,k,l
 *
 * hkl = UBI.gvector
 */
int hklcalc( double XL[], double axis[], double ang,
      double UBI[], double T[],
      double dLn[], double O[], double *M,
      double kvector[],
      double gvector[], 
      double hkl[],
      double wvln ){
    int err;

    err = rotate_vector_axis_angle( axis, -ang, T, O ); 
    if ( err != 0 ) return err;

    /* Unit vector in lab is: */
    dLn[0] = XL[0] - O[0];
    dLn[1] = XL[1] - O[1];
    dLn[2] = XL[2] - O[2];
    *M = norm3(dLn);
    err = normalise_vector( dLn );
    if ( err != 0 ) return err;
    
    kvector[0] = (dLn[0] - 1.0)/wvln; /* Beam along X */
    kvector[1] = (dLn[1] )/wvln;
    kvector[2] = (dLn[2] )/wvln;

    err = rotate_vector_axis_angle( axis, ang, kvector, gvector );
    if ( err != 0 ) return err;
    
    mat3_transform_vec( UBI, gvector, hkl );

    return 0;
}

/**
 * Determine the omega angle which would minimise the error
 *
 * Perhaps... numerical derivative followed by dot product on error in h?
 * Or ... solve the inverse equations to get omega?
 */

/**
 * Compute the apparent grain origin backwards from the integer
 * h,k,l indices
 *
 * Assume the distance of detector to grain is preserved, which gives
 * the effective 2theta error part
 */
int computeT( double XL[], double ih[], double UB[], double gcalc[], double axis[],
        double ang, double kcalc[], double M, double dLcalc[], double Ocalc[],
        double Tcalc[], double wvln){
    int err;

    mat3_transform_vec( UB, ih, gcalc);
    
    err = rotate_vector_axis_angle( axis, -ang, gcalc, kcalc );

    if ( err !=0 ) return err;

    dLcalc[0] = (kcalc[0] * wvln + 1) * M;
    dLcalc[1] = kcalc[1] * wvln * M;
    dLcalc[2] = kcalc[2] * wvln * M;

    Ocalc[0] = XL[0] - dLcalc[0];
    Ocalc[1] = XL[1] - dLcalc[1];
    Ocalc[2] = XL[2] - dLcalc[2];
    
    err = rotate_vector_axis_angle( axis, ang, Ocalc, Tcalc );
    return err;
}




void printvec3( char* name, double vec[] ){
    printf("%s %f %f %f\n",name,vec[0],vec[1],vec[2]);
}

int main(){
    /* inputs */
    double XL[3]={100.,100.,100.};
    double axis[] = {0.,0.,1.};
    double ang = 0.3;
    double UBI[9] = { 1.,0.,0.,0.,1.,0.,0.,0.,1.};
    double T[3] = {10,11,12};
    double wvln = 0.154;
    /* outputs */
    double kvector[3], gvector[3], hkl[3], dL[3], O[3], M, UB[9];
    double kcalc[3], gcalc[3],dLcalc[3],Ocalc[3],Tcalc[3];
    double ih[3];
    int err;

    printf("Hello\n");

    hklcalc(XL, axis, ang, UBI, T, dL, O, &M, kvector, gvector, hkl, wvln );
    printvec3("dL",dL);
    printf("M %f\n",M);
    printf("hkl %f %f %f ",hkl[0],hkl[1],hkl[2]);
    round3( hkl, ih );
    printf("hkl %f %f %f \n",ih[0],ih[1],ih[2]);
    
    err = inverse_mat3( UBI, UB );
    if (err != 0) return err;

    computeT(XL, ih, UB, gcalc, axis, ang, kcalc, M, dLcalc, Ocalc, Tcalc, wvln );

    printvec3("Gcalc",gcalc);
    printvec3("Kcalc",kcalc);
    printvec3("dLcalc",dLcalc);
    printf("M %f\n",M);
    printvec3("Tcalc",Tcalc);


    return 0;
}
