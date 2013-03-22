
#include <math.h> /* fabs */
#include "math.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif
/* Affine transformations */


/* Functions 
 * Return value is error status
 */

typedef double real;
typedef real rmatrix[9];
typedef real vector[3];
typedef struct {
    rmatrix m;
    vector  v;
} afmatrix;

typedef real scalar;


/* Written and at least 1 test case */
void mat3_transform_vec( rmatrix, vector, vector*);
void mat3_T_transform_vec( rmatrix m, vector v, vector *r);

scalar determinant3( rmatrix );

int inverse_mat3( rmatrix, rmatrix*);

void affine_transform_vec( afmatrix, vector, vector*);

void mat3_prod( rmatrix , rmatrix, rmatrix*);

int inverse_affine_mat( afmatrix, afmatrix*);

int normalise_vector(vector*);

scalar norm3(vector);

int mat3_from_axis_angle(vector*, scalar , rmatrix* );

/* Written and untested */
double degrees(double);
double radians(double);

#define _norm3( a ) \
    (sqrt( (a)[0]*(a)[0] + (a)[1]*(a)[1] + (a)[2]*(a)[2] ))

#define _dot3( a, b )\
    (( (a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2] ))

#define round3( a, b ) \
    (b)[0] = floor( a[0] + 0.5); \
    (b)[1] = floor( a[1] + 0.5); \
    (b)[2] = floor( a[2] + 0.5); 


#define crossProduct(a,b,c) \
	(c)[0] = (a)[1] * (b)[2] - (b)[1] * (a)[2]; \
	(c)[1] = (a)[2] * (b)[0] - (b)[2] * (a)[0]; \
	(c)[2] = (a)[0] * (b)[1] - (b)[0] * (a)[1];

#define vec3sub(a,b,c) \
	(c)[0] = (a)[0] - (b)[0]; \
	(c)[1] = (a)[1] - (b)[1]; \
	(c)[2] = (a)[2] - (b)[2]; \

#define vec3sdiv(a,b,c) \
    (c)[0] = (a)[0]/(b);\
    (c)[1] = (a)[1]/(b);\
    (c)[2] = (a)[2]/(b);\



#define vec3add(a,b,c) \
	(c)[0] = (a)[0] + (b)[0]; \
	(c)[1] = (a)[1] + (b)[1]; \
	(c)[2] = (a)[2] + (b)[2]; \

#define matVec(a,b,c) \
    (c)[0] = (a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2]; \
    (c)[1] = (a)[3]*(b)[0] + (a)[4]*(b)[1] + (a)[5]*(b)[2]; \
    (c)[2] = (a)[6]*(b)[0] + (a)[7]*(b)[1] + (a)[8]*(b)[2]; 

#define matTVec(a,b,c) \
    (c)[0] = (a)[0]*(b)[0] + (a)[3]*(b)[1] + (a)[6]*(b)[2]; \
    (c)[1] = (a)[1]*(b)[0] + (a)[4]*(b)[1] + (a)[7]*(b)[2]; \
    (c)[2] = (a)[2]*(b)[0] + (a)[5]*(b)[1] + (a)[8]*(b)[2]; 

int rotate_vector_axis_angle( vector*, scalar , vector, vector *);

/* Yet to be written */





/* Derivatives of the previous functions ! */
double d_determinant3( double [], double []);
/* Matrix inverse 
 * There are 9 outputs and 9 inputs
 * Therefore 81 derivatives d_out/d_in
 * */

int omegacalc( vector ,  rmatrix , rmatrix ,
        vector , real , real* , real* );

int hklcalc( vector , vector *, real ,
      rmatrix , vector ,
      vector *, vector *, real *,
      vector *,
      vector *, 
      vector *,
      real );

void printvec3( char* , real* );









