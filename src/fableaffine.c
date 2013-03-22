#include "affine.h"
#include <stdio.h>





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
