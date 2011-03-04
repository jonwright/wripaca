/* wripaca.h
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

#define REAL double

/**
 * Linear regular binning functions
 * origin, step
 * For log binning take log(Q) from two theta etc
 */
typedef struct {
    REAL origin; /* Point at i==0 */
    REAL step;
    int npt;
    int circular;
} bins;

int getbin(bins *, REAL );

REAL getboundary(bins *, int );


/**
 * A polygon datatype. List of points, always in order to walk 
 * around edges
 *
 * src = The source image co-ordinates. i,j or spatial positions
 * des = The destination image co-ordinates. tth eta etc
 * bin = The bin which a pixel belongs to
 * maxb = The maximum bin number the polygon touches (for looping)
 * np = The number of points making the polygon, max (NP==8)
 * first = The smallest bin number point in the polygon
 */
#define NP 8
#define NDIM 2

typedef struct {
    REAL src[NP*NDIM]; /* Point pos in src image */
    REAL des[NP*NDIM]; /* Point pos in dest image */
    int bin[NP*NDIM];    /* Point bin in dest image */
    int maxb[NDIM];      /* Highest bin values cache */
    int np;           /* Number of points */
    int first;        /* Point to start looping from */
} polygon;


REAL polygon_area( polygon * );

void polygon_compute_corners( polygon *, 
        REAL (*)(REAL []), 
        REAL (*)(REAL [])
        );

void polygon_compute_bins( polygon *, 
       bins *,
       int  );

void polygon_setfirst( polygon *, int );

int polygon_cut_edge( polygon *, int , int , int , REAL , 
        polygon *);

int polygon_cut( polygon *, polygon *, polygon *,
        int , int , int , REAL );


void polygon_create(polygon *, int , int ,
        REAL (*)(REAL []),
        REAL (*)(REAL []));


void test();


