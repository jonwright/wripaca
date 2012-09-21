

#include "affine.h"
#include "math.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

/**
 * Convert degrees to radians
 *
 * @param angle in degrees
 * @return angle in radians
 */
double radians(double p){ return p*M_PI/180.0; }

/**
 * Convert radians to degrees
 *
 * @param angle in radians
 * @return angle in degrees
 */
double degrees(double p){ return p*180.0/M_PI; }


/**
 * Round to nearest integers
 *
 * @param double vector
 * @return nearest integers
 */

void round3(double v[], double h[]){
    h[0] = floor(v[0]+0.5);
    h[1] = floor(v[1]+0.5);
    h[2] = floor(v[2]+0.5);
}

/**
 * Vector 3 dot product
 *
 * @param a
 * @param b
 * @return a.b
 */
scalar dot3( vector a, vector b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/**
 * Vector 3 cross product
 *
 * @param a
 * @param b
 * @return a x b
 */
void cross3(vector a, vector b, vector *c){
    (*c)[0] = a[1]*b[2] - a[2]*b[1];
    (*c)[1] = a[2]*b[0] - a[0]*b[2];
    (*c)[2] = a[0]*b[1] - a[1]*b[0];
}
/**
 * Computes a 3x3 matrix vector product 
 *
 * @param double[9] m is a 3x3 matrix stored flat in row order
 * @param double[3] v is a vector
 * @param double[3] r is the result
 */
void mat3_transform_vec( rmatrix m, vector v, vector *r){
    (*r)[0] = m[0] * v[0] + m[1]*v[1] + m[2]*v[2];
    (*r)[1] = m[3] * v[0] + m[4]*v[1] + m[5]*v[2];
    (*r)[2] = m[6] * v[0] + m[7]*v[1] + m[8]*v[2];
}


/**
 * Computes 3x3 matrix product c = a.b
 *
 * All matrices are 3x3 stored in row major order
 * 0 1 2  a b c   0a 1d 2g, 0b 1e 2h,  etc
 * 3 4 5  d e f = 
 * 6 7 8  g h i   
 *
 * @param double[9] a 
 * @param double[9] b
 * @param double[9] c 
 * @return error status
 */
void mat3_prod( rmatrix a, rmatrix b, rmatrix *c){
    (*c)[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
    (*c)[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
    (*c)[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
    (*c)[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
    (*c)[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
    (*c)[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
    (*c)[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
    (*c)[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
    (*c)[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}

/**
 * Does an affine transform matrix vector product 
 *
 * The affine matrix + vector is supplied in a struct
 * as (m,v)
 *
 * @param double[12] m is a 4x3 matrix stored in row order
 * @param double[3] v is a vector
 * @param double[3] r is the result
 * @return error status 
 */
void affine_transform_vec( afmatrix m, vector v, vector *r){
    (*r)[0] = m.m[0]*v[0] + m.m[1]*v[1] + m.m[2]*v[2] + m.v[0];
    (*r)[1] = m.m[3]*v[0] + m.m[4]*v[1] + m.m[5]*v[2] + m.v[1];
    (*r)[2] = m.m[6]*v[0] + m.m[7]*v[1] + m.m[8]*v[2] + m.v[2];
}


/**
 * Finds the determinant of a 3x3 matrix
 *
 * @param rmatrix is a 3x3 matrix with row major storage
 * @return determinant
 */
scalar determinant3( rmatrix m ){
    return m[0]*m[4]*m[8] - m[0]*m[5]*m[7] + m[1]*m[5]*m[6] - 
           m[1]*m[3]*m[8] + m[2]*m[3]*m[7] - m[2]*m[4]*m[6];
}


/**
 * Inverts a 3x3 matrix
 *
 * |0 1 2 |-1             | 
 * |3 4 5 |    =  1/DET * | 
 * |6 7 8 |               | <
 *
 * @param m is a 3x3 matrix with row major storage
 * @param r will be the inverse 
 * @return error if the determinant is zero
 */
int inverse_mat3(rmatrix m, rmatrix *r){
    scalar d;
    d = determinant3( m );
    if (fabs(d) == 0.0F){  
        return 1;
    } 
    d = 1.0f/d;
    (*r)[0] =  m[4]*m[8]*d - m[7]*m[5]*d;
    (*r)[1] =  m[7]*m[2]*d - m[8]*m[1]*d;
    (*r)[2] =  m[1]*m[5]*d - m[2]*m[4]*d;
    (*r)[3] =  m[5]*m[6]*d - m[3]*m[8]*d;
    (*r)[4] =  m[8]*m[0]*d - m[2]*m[6]*d;
    (*r)[5] =  m[2]*m[3]*d - m[0]*m[5]*d;
    (*r)[6] =  m[3]*m[7]*d - m[6]*m[4]*d;
    (*r)[7] =  m[1]*m[6]*d - m[0]*m[7]*d;
    (*r)[8] =  m[0]*m[4]*d - m[3]*m[1]*d;
    return 0;
}


/**
 * Find the inverse of an affine transform matrix
 *
 * r = A.x + b
 * x = A-1.(r - b)
 * x = A-1.r - A-1.b
 *
 * @param m is the affine transform matrix 12= 3x3 + 3 trans
 * @param r is the resulting inverse
 * @return error if determinant is zero
 */
int inverse_affine_mat(afmatrix m, afmatrix *r){
    int err;
    err = inverse_mat3( m.m, &(*r).m );
    if(err != 0){
        return 1;
    }
    mat3_transform_vec((*r).m, m.v, &(*r).v);
    (*r).v[0] = -(*r).v[0];
    (*r).v[1] = -(*r).v[1];
    (*r).v[2] = -(*r).v[2];
    return 0;
}

/**
 * Find the length of a vector
 *
 * @param a length 3 vector
 * @return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
 */
scalar norm3(rmatrix v){
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

/**
 * Normalise a vector to unit length
 *
 * @param a length 3 double vector
 * @return error if the vector is zero
 */
int normalise_vector(vector *v){
    scalar norm;
    norm = norm3(*v);
    if( norm > 0 ){
       (*v)[0] = (*v)[0]/norm;
       (*v)[1] = (*v)[1]/norm;
       (*v)[2] = (*v)[2]/norm;
       return 0;
    }
    return 1;
}

/**
 * Find the 3x3 rotation matrix from axis and angle description
 *
 * [R] =	0 t*x*x + c	, 1 t*x*y - z*s	, 2 t*x*z + y*s
 *          3 t*x*y + z*s, 4 t*y*y + c	    , 5 t*y*z - x*s
 *          6 t*x*z - y*s	, 7 t*y*z + x*s	, t*z*z + c
 *
 * http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToMatrix/index.htm
 *
 * @param a = axis - must be able to normalise this
 * @param t = angle in radians
 * @param m = the matrix to be computed
 * @return error if the vector cannot be normalised
 */
int mat3_from_axis_angle(vector *a, scalar p, rmatrix *m){
    int err;
    double c,s,t,T1,T2;
    err = normalise_vector(a);
    if (err != 0){ return err;}
    s = sin(p);
    c = cos(p);
    t = 1-c;
    (*m)[0] = t*(*a)[0]*(*a)[0] + c;
    (*m)[4] = t*(*a)[1]*(*a)[1] + c;
    (*m)[8] = t*(*a)[2]*(*a)[2] + c;
    T1 = t*(*a)[0]*(*a)[1];
    T2 = (*a)[2]*s;
    (*m)[1] =  T1 - T2;
    (*m)[3] =  T1 + T2;
    T1 = t*(*a)[1]*(*a)[2];
    T2 =  (*a)[0]*s;  
    (*m)[5] = T1 - T2;
    (*m)[7] = T1 + T2;
    T1 = t*(*a)[0]*(*a)[2]; 
    T2 = (*a)[1]*s;
    (*m)[2] = T1 + T2;
    (*m)[6] = T1 - T2;
    return 0;
}


/**
 * Apply the axis/angle formula to rotate a vector an axis
 *
 * See  http://mathworld.wolfram.com/RotationFormula.html
 * r' = r cos(p) + n.(n.r)(1-cos(p)) + (rxn).sin(p)
 * Implemented mainly for testing the matrix formulation 
 *
 * @param a[3] is the axis direction
 * @param p is the angle in radians
 * @param v[3] is the input point
 * @param r[3] is the ouput point
 */
int rotate_vector_axis_angle(vector *a, scalar p, vector v, vector *r){
    scalar sinp, cosp, tmp;
    vector t;
    int err;
    err = normalise_vector(a);
    if (err != 0){ return err;}
    cosp = cos(p);
    sinp = sin(p);
    tmp = dot3((*a),v)*(1-cosp);
    cross3(v,(*a),&t);
    (*r)[0] = v[0] * cosp + tmp*(*a)[0] + t[0]*sinp ;
    (*r)[1] = v[1] * cosp + tmp*(*a)[1] + t[1]*sinp ;
    (*r)[2] = v[2] * cosp + tmp*(*a)[2] + t[2]*sinp ;
    return 0;
}


