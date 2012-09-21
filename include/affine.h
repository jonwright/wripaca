
#include <math.h> /* fabs */
/* Affine transformations */


/* Functions 
 * Return value is error status
 */

typedef float rmatrix[9];
typedef float vector[3];
typedef struct {
    rmatrix m;
    vector  v;
} afmatrix;

typedef float scalar;


/* Written and at least 1 test case */
void mat3_transform_vec( rmatrix, vector, vector*);

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
void round3( double[], double[]);
scalar dot3( vector, vector);
void cross3( vector, vector, vector *);

int rotate_vector_axis_angle( vector*, scalar , vector, vector *);

/* Yet to be written */





/* Derivatives of the previous functions ! */
double d_determinant3( double [], double []);
/* Matrix inverse 
 * There are 9 outputs and 9 inputs
 * Therefore 81 derivatives d_out/d_in
 * */














