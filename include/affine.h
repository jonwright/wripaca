
#include <math.h> /* fabs */
/* Affine transformations */


/* Functions 
 * Return value is error status
 */

/* Written and at least 1 test case */
void mat3_transform_vec( double [], double [], double []);

double determinant3( double[]);

int inverse_mat3( double [], double []);

void affine_transform_vec( double [], double [], double []);

void mat3_prod( double [], double[], double[]);

int inverse_affine_mat( double [], double []);

int normalise_vector(double []);

double norm3(double []);

int mat3_from_axis_angle(double [], double , double []);

/* Written and untested */
double degrees(double);
double radians(double);
void round3( double[], double[]);
double dot3( double[], double[]);

int rotate_vector_axis_angle( double [], double , double [], double []);

/* Yet to be written */





/* Derivatives of the previous functions ! */
double d_determinant3( double [], double []);
/* Matrix inverse 
 * There are 9 outputs and 9 inputs
 * Therefore 81 derivatives d_out/d_in
 * */














