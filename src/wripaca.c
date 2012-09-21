/* wripaca.c
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


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "wripaca.h"


/* FIXME
 * eta cut points
 * beam centre containing pixel
 * binning outside of limits - clipping bins == -1
 * binning as a function of the other dimensions bin (Peter's fazit idea)
 */


/**
 * To deal with problems in a cleaner way for gui etc
 * This prints when the error happened, like assert does
 *
 * (strongly) inspired by P. J. Plauger's book The C standard library.
 */
#define _STR(x) _VAL(x)
#define _VAL(x) #x
#define IFOK( test )  ( (~test) ? (void) 0 \
   : handle_error(__FILE__ ":" _STR(__LINE__) " " #test) )

/**
 * Private function definitions
 */
int polygon_append(polygon *, polygon *, int );
int polygon_iter( polygon *, int );
/* The library shall not print */
void printpoint( polygon *, int , FILE*  );
void polygon_print( polygon *, FILE*, char * );


void polygon_plot(polygon *,  int , int , FILE* );

/* Marginally faster 
#define IFOK( test ) test 
*/

void handle_error(char *msg){
    fprintf(stderr,"Error! Giving up, sorry\n");
    fputs(msg, stderr);
    exit(1);
}


/**
 * Finds out which bin a REAL value is in
 * Uses -1 for all bins which are out of range (to be caught later)
 *
 * @param bins - binning structure (origin, step, npts)
 * @param v the value to bin
 * @returns the bin value
 */
int getbin(bins *b, REAL v){
    int ib;
    ib = (int) floor( (v - b->origin)/b->step) ;
    if(b->circular) return ib % (b->npt); /* Modulo nsteps */
    return ib;
}


/**
 * Finds the lower boundary of bin i
 *
 * @param bins - binning structure (origin, steps, npts, circular)
 * @param i - which bin
 * @returns the lower boundary of bin i
 */
REAL getboundary(bins *b, int i){
    return b->origin + b->step * i;
}

/**
 * Finds the direction to iterate
 * FIXMi
 */


/**
 * Appends a point onto a polygon
 * 
 * @param polygon dest to append point onto
 * @param polygon src to read point from
 * @param i is the point in src to append
 * @returns 0 for success, 1 for failure (due too many points, 
 * should not happen in normal usage)
 */
int polygon_append(polygon *dest, polygon *src, int i){
    int j,k;
    k = dest->np; /* Last point in dest */
    dest->np++; /* Add new point */
    if( dest->np >= NP ){
        return 1; /* Too many points */
    }
    for(j=0; j<NDIM; j++){ /* Copy i in src to k in dest */
       dest->src[NDIM*k+j] = src->src[i+j];
       dest->des[NDIM*k+j] = src->des[i+j];
       dest->bin[NDIM*k+j] = src->bin[i+j];
    }
    return 0;
}



/**
 * Helper function to iterate around a polygon
 * Intended to be inlined
 *
 * @param polygon p which is being iterated over
 * @param i is the current position on the polygon
 */
int polygon_iter( polygon *p, int i){
    assert(i%NDIM == 0);
    assert(i <= NDIM*p->np);
    i += NDIM;
    if ( i == NDIM*p->np ) i = 0;
    return i;
}

/** 
 * Prints a single point from a polygon to stdout
 *
 * @param polygon p holding the point
 * @param i is which point to print
 */
void printpoint( polygon *p, int i, FILE* output ){
    int j;
    printf("%2d ",i);
    assert( i % NDIM == 0 );
    assert( (i >= 0) && (i <= NDIM*p->np));
    fprintf(output,"%2d ",i/NDIM);
    for(j=0;j<NDIM;j++){
        fprintf(output,"s%d %8.4f d%d %8.4f b%d %4d  ", 
                        j,p->src[i+j],
                        j,p->des[i+j],
                        j,p->bin[i+j]);
    }
    fprintf(output,"\n");
}


/**
 * Prints a polygon to stream 
 *
 * @param  polygon p to be printed
 * @param  output is a stream to print to (stdout if NULL)
 */
void polygon_print( polygon *p, FILE* output, char * name ){
    int i;
    if(output == NULL) {
        return;
    }
    assert(p->first % NDIM == 0);
    fprintf(output,"Polygon %s, npoints = %d, first = %d\n",name,
            p->np,p->first/NDIM);
    if( p-> np > 0 ){
        i = p->first;
        do{
            printpoint(p, i, output);
            i = polygon_iter( p, i );
        } while ( i != p->first );
    }
}

/**
 * Computes the area of a polygon by walking around the edges
 *
 * @param polygon p
 */
REAL polygon_area( polygon *p ){
    /* FIXME - how does this work for NDIM != 2 ?? */
    REAL a, x0, y0;
    int last, next;
    assert(NDIM == 2);
    a=0;
    last = p->first;
    x0 = p->src[p->first]; /* We keep these to subtract for stability */
    y0 = p->src[p->first+1];
    do{
        next = last+2;
        if (next == 2*p->np) next = 0;
        a += (p->src[last]-x0) * (p->src[next+1] - y0);
        a -= (p->src[next]-x0) * (p->src[last+1] - y0);
        last = next;
    } while ( last != p->first );
    return a/2.0;
}

/**
 * Obtains the corner destination values for the current polygon
 * We normally plan to call this once per pixel for the four corners only
 *
 * @param p polygon with points defined
 * @param get0 is the function for binning dimension 0
 * @param get1 is the function for binning dimension 1
 */
void polygon_compute_corners( polygon *p, 
        REAL (*get0)(REAL []), 
        REAL (*get1)(REAL [])
        ){
    int i;
    i = p->first;
    do{
        p->des[i]   = (*get0)(&p->src[i]);
        p->des[i+1] = (*get1)(&p->src[i]);
        i = polygon_iter( p, i );
    } while ( i != p->first );
}

/**
 * Decides which bin a point is in, given the values in p->des
 *
 * @param polygon p to do computation on
 * @param binning structure to use for bins
 * @param dim the dimension which is being binned
 */
void polygon_compute_bins( polygon *p, 
       bins *b,
       int dim ){
    int i;
    assert(dim >= 0);
    assert(dim < NDIM);
    i = p->first;
    do{
        p->bin[i+dim]   = getbin(b, p->des[i+dim]);
        i = polygon_iter( p, i );
    } while ( i != p->first );
}

/**
 * Decide which point to start interating around the polygon
 * The first point should correspond to the first (lowest) bin
 * Fills in p->first and also p->maxb[dim]
 *
 * @param p is the polygon to look at
 * @param dim is the dimension which is being searched
 */
void polygon_setfirst( polygon *p, int dim){
    int i, minb, mini, maxb, maxi;
    /* Initialise with the current first bin */
    i = p->first;
    minb = p->bin[i+dim];
    mini = i;
    maxb = p->bin[i+dim];
    maxi = i;
    do{
        if( p->bin[i+dim] < minb ){
            minb =  p->bin[i+dim];
            mini = i;
        }
        if( p->bin[i+dim] > maxb ){
            maxb =  p->bin[i+dim];
            maxi = i;
        }
        i = polygon_iter( p, i );
    } while ( i != p->first );
    p->first = mini;
    p->maxb[dim] = maxb;
    /* Out of range is mini == -1 */
}

/**
 * Cuts edge a -> b in polygon p at bound and writes to end of dest
 *
 * @param src is the polygon containing an edge to be cut
 * @param a is the start point for cutting
 * @param b is the end point for cutting
 * @param dim is the dimension which is being cut (eg: either tth or eta)
 * @param bound is the value at which to cut
 * @param dest is the polygon where the new point is addeded
 */
int polygon_cut_edge( polygon *src, int a, int b, int dim, REAL bound, 
        polygon *dest){
    REAL wa, wb, diff;
    int j,k;
    /* FIXME : check bound is between a and b */
    /* The length of the edge being cut */
    diff = src->des[a+dim] - src->des[b+dim];
    if(diff == 0){
        /* FIXME This might happen, what to do? */
        printf("Error in polygon_cut_edge: zero length edge %d %d\n",a,b);
        polygon_print(src, stderr,"src in cut_edge");
        exit(0);
    }
    wa = (bound - src->des[b+dim])/diff;
    wb = (src->des[a+dim] - bound)/diff;
    /* Append on dest */
    k = dest->np;
    dest->np++;
    if( k >= NP ) 
        return 1; 
    for(j=0; j<2; j++){
        dest->src[2*k+j] = wa * src->src[a+j] + wb * src->src[b+j];
        dest->des[2*k+j] = wa * src->des[a+j] + wb * src->des[b+j];
    }
    return 0;
}



/**
 * Cuts the 'thisbin' corner off of polygon p putting the ear in e and the 
 * remainded in r
 *
 * @param the polygon to have an ear cut off
 * @param e is the returned ear (overwritten)
 * @param r is the remainder after ear removal (overwritten)
 * @param is the dimension being cut, eg tth or eta
 * @param thisbin is the target bin for the ear
 * @param nextbin is the next bin for the cutline
 * @param bound is the value at which to cut
 */
int polygon_cut( polygon *p, polygon *e, polygon *r,
        int dim, int thisbin, int nextbin, REAL bound){
    int inlo, i, prev, iend;
    /* We start out inside the ear */
    inlo = 1;
    i = p->first;
    assert( p->bin[p->first + dim] == thisbin); /* We are in the low bin */
    /* Initialise the output polygons for overwriting */
    e->np = 0; e->first = 0;
    r->np = 0; r->first = 0;
    /* This point goes straight into the ear */
    IFOK( polygon_append(e, p, i) );
    prev = i; /* prev -> i is the current edge */
    i = polygon_iter(p, i);
    iend = i; /* Stop the loop here */
    do {
        if( inlo ){ /* In the ear being cut off */
            if( p->bin[i+dim] == thisbin ){
                if ( i != p->first ){ /* Dont add when re-processing at end */
                    IFOK( polygon_append(e, p, i) );
                }
            } else { /* Switch from thisbin to nextbin here */
                /* Appends new point at end of e */
                IFOK( polygon_cut_edge( p, prev, i, dim, bound, e ));
                /* Now, e has a new point, should be edge of thisbin */
                e->bin[2*(e->np-1) + dim] = thisbin;
                /* Copy the new point to the remainder and set the bin */
                IFOK( polygon_append(r, e, 2*(e->np-1)));
                r->bin[2*(r->np-1) + dim] = nextbin;
                if( i != p->first) {
                    /* Add also the other end of the edge, if we are not back
                     * at the start */
                    IFOK (polygon_append( r, p, i ));
                }
                inlo = 0; /* Now we are walking along the high edge */
            }
        } else { /* not inlo, now we are inhi */
            if( p->bin[i+dim] != thisbin ){
                IFOK( polygon_append(r, p, i));
            } else { /* Switch back from hibin to thisbin */\
                IFOK( polygon_cut_edge(p, prev, i, dim, bound, e));
                /* Appends the new point directly on e (joins last appended)*/
                e->bin[2*(e->np-1) + dim] = thisbin;
                /* And copy it to r too */
                IFOK(polygon_append(r, e, 2*(e->np-1)));
                r->bin[2*(r->np-1) + dim] = nextbin;
                if( i != p->first ) { 
                    /* Append the end of this edge if not already done */
                    IFOK(polygon_append( e, p, i ));
                }
                inlo = 1;
            }
        }        
        prev = i;
        i = polygon_iter(p, i);
    } while ( i != iend );
    /* Here we check if we made a single cut or not */
    if ((e->np + r->np) == 2+(p->np))  return 0;
    else return 1;
}

/**
 * Plots a polygon at pixel i,j - used now for printing and debug
 *
 */
void polygon_plot(polygon *p,  int i, int j, FILE* outfile ){
    int k;
    assert(p->first % 2 == 0);
    if( p-> np > 0 ){
        fprintf(outfile,"\n\n");
        k = p->first;
        do{
            fprintf(outfile,
                    "%8.4f %8.4f %d %d \n",p->src[k],p->src[k+1],i,j);
            k = polygon_iter( p, k );
        } while ( k != p->first );
    }
}

/**
 * Initialise a polygon for test case
 *
 * When using the code for a spatial distortion this function
 * would look up the spatial values at the corners
 * FIXME: This should be an argument to dopixel
 *
 * @param polygon to be created (overwritten)
 * @param i is the slow pixel co-ordinate
 * @param j is the fast pixel co-ordinate
 * @param get0 is a function to get the two theta (etc) value
 * @param get1 is a function to get the azimuth value
 */
void polygon_create(polygon *o, int i, int j,
        REAL (*get0)(REAL []),
        REAL (*get1)(REAL [])){
    /* Create the new polygon */
    o->src[0] = -0.5 + i ; o->src[1]=-0.5 + j;
    o->src[2] =  0.5 + i ; o->src[3]=-0.5 + j;
    o->src[4] =  0.5 + i ; o->src[5]= 0.5 + j;
    o->src[6] = -0.5 + i ; o->src[7]= 0.5 + j;
    o->first = 0;
    o->np = 4;
    /* printf("Area = %f\n",polygon_area(&o)); */
    /* Initialise the corner data */
    polygon_compute_corners( o, get0, get1);
}



/**
 * Process a single pixel i,j 
 *
 * @param i is the i co-ordinate
 * @param j is the j co-ordinate
 * @param get0 is a function to compute, eg, tth
 * @param get1 is a function to compute, eg, eta
 * @param b0 is the binning structure for get0
 * @param b1 is the binning structure for get1
 * @outfile is an outputfile is you are printing debug info
 */
int dopixel( int i, int j,
        REAL (*get0)(REAL []),
        REAL (*get1)(REAL []),
        bins *b0,
        bins *b1){

    polygon *o,*e0,*r0,*e1,*r1,*t; /* Ear-0, Rest-0,
                              * Ear-1, Rest-1 */
    polygon oa,e0a,r0a,e1a,r1a; /* ..a allocated */
    int i1, i0;
    int  max0, max1, np;
    REAL val, areasum;

    /* Initialise the destination polygons and sums */
    e0a.np = 0; r0a.np = 0; e1a.np =0 ; r1a.np =0;
    e0 = &e0a ; e1 = &e1a; r0 = &r0a; r1=&r1a; o = &oa;
    areasum = 0;
    np = 0;

    /* Initialise the source pixel */
    polygon_create(o, i, j, get0, get1);
    /* Decide on bins for looping dim0 */
    polygon_compute_bins( o, b0, 0);
    /* Check if all points are out of bounds */
    polygon_setfirst( o, 0 );
    max0 = o->maxb[0];
    /* Out of bounds for all points in polygon */
    if( max0 < 0) return 0;
    if( o->bin[o->first] > b0->npt ) return 0;
    /* At least something in bounds */
    for( i0 = o->bin[o->first]; i0 <= max0; i0++){ /* Loop first dim */
        val = getboundary( b0, i0+1 );
        /* FIXME nextbin for eta */
        IFOK( polygon_cut( o, e0, r0, 0, i0, i0+1, val ) );
        /* Now decide looping for dim1 */
        polygon_compute_bins( e0, b1 , 1);
        polygon_setfirst( e0, 1 );
        max1 = e0->maxb[1];
        if( max1 < 0 ) continue; /* Out of bounds */
        if( e0->bin[e0->first+1] > b1->npt ) continue; /* Out of bounds */
        /* FIXME: Out of azimuth or two theta bin range ?? */
        for( i1 = e0->bin[e0->first+1]; i1<=max1; i1++){
            val = getboundary( b1, i1+1); /* FIXME looparound at 180 -> -180 */
            IFOK (polygon_cut( e0, e1, r1, 1, i1, i1+1, val) );
            /* OUTPUT e1 */
            areasum += polygon_area(e1);
            np++;
			t = e0; /* swap */
			e0 = r1;
			r1 = t;
        }
		t = o; /* swap */
		o = r0;
		r0 = t;
    } /* Loop i0 - the tth bins */
    return np;
}


