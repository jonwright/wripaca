/* testwripaca.c
 *
 *
 * Copyright (C) 2011  Jon Wright
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  0211-1307  USA
 */


#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "wripaca.h"

/**
 * Model tth function as stand-in for testing
 * Distance is 100 mm, center is 0,0
 *
 * @param x[2] is the x/y position on a detector
 * @returns two theta in degrees
 */
REAL gettth(REAL x[]){
   REAL val;
   /* printf("tth x y %f %f\n",x[0],x[1]); */
   val = atan2( sqrt(x[0]*x[0]+x[1]*x[1]), 100.0) * 180.0/3.14159;
   return val;
}

/**
 * Model eta function as stand in for testing
 * center is 0,0
 *
 * @param x[2] is the x/y position on detector
 * @returns eta in degrees
 */
REAL geteta(REAL x[]){
    return atan2( x[0], x[1] ) * 180.0/3.14159;
}

/**
 * Uses openmp to find the time to process a 2K by 2K image
 */
int testforspeed( bins *bs, bins *bf, 
        REAL (*gettth)(REAL[]), 
        REAL (*geteta)(REAL[])){
    int np, nt, nthreads, i,j;
    clock_t start, end;
    REAL elapsed;
    np = 2048;
    nt = 0;
    nthreads = 1;
    start = clock();
#pragma omp parallel for private(j) reduction(+:nt)
    for(i=1;i<np;i++){
#ifdef _OPENMP
    if( i==1 ) nthreads = omp_get_num_threads();
#endif
        for(j=1;j<np;j++){
            nt = nt + dopixel(i,j, gettth, geteta, bs, bf);
        }
    }
    end = clock();
    elapsed = ((REAL) end - (REAL)start)/((REAL) CLOCKS_PER_SEC);
    printf("Cutting pixels took %f, %f ms per pixel\n",elapsed,1e3*elapsed/np/np);
    printf("Got %d fragments, %f per pixel\n",nt,((REAL)nt)/np/np);
    printf("Used %d threads\n",nthreads);
    return 0;
}

/**
 * For a specific pixel i, j this prints the segments to files
 * for plotting with gnuplot
 */
int testforgnuplot( bins *bs, bins *bf, 
        REAL (*gettth)(REAL[]), 
        REAL (*geteta)(REAL[]), int i, int j){
    int n;
    clock_t start, end;
    REAL elapsed;
    static FILE* outfile;
    printf("Example will fail!");
    exit(1);
    outfile = fopen( "surfbin.dat", "w" );
    start = clock();
    n = dopixel(i,j, gettth, geteta, bs, bf);
    end = clock();
    elapsed = ((REAL)end - (REAL)start)/(REAL) CLOCKS_PER_SEC;
    fclose(outfile);
    outfile = fopen( "gnuplot.input","w");
    fprintf(outfile,"set term win\nset nokey\n");
    fprintf(outfile,"plot \"surfbin.dat\" i %d u 1:2 w filledcurves closed\n",0);
    for(i=1;i<n;i++)
        fprintf(outfile,
                "replot \"surfbin.dat\" i %d u 1:2 w filledcurves closed\n",
                i);
    fclose(outfile);
    return 0;
}

void test(){
    bins bs = { 0, 0.1, 90000, 0};
    bins bf = { -360,  1, 720, 1 };

    testforspeed( &bs, &bf, gettth, geteta);
}

int main(int argc, char* argv[]){
        test();
        return 0;
}

/* http://www.azillionmonkeys.com/qed/hash.html */
