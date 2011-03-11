
/* Unit test code from minunit
 * http://pantras.free.fr/articles/minunit.html
 * ===========================================================================
 */

#define MU_ASSERT(message, test) do { if (!(test)) return message; } while (0)
#define MU_RUN_TEST(test) do { char *message = test(); tests_run++; \
                               if (message) return message; } while (0)
int tests_run=0;

#include "affine.h"
#include <stdio.h>




/* Test functions follow
 * ===========================================================================
 */


char* test_mat3_transform_vec(void){
    double vec[3] = { 1, 2, 3 };
    double mat[9] = { 0, 0, 1, 
                     0, 1, 0, 
                     1, 0, 0 };
    double t[3];
    mat3_transform_vec( mat, vec, t );
    MU_ASSERT( "transform_vec0", t[0] == vec[2] );
    MU_ASSERT( "transform_vec1", t[1] == vec[1] );
    MU_ASSERT( "transform_vec2", t[2] == vec[0] );
    return 0;
 
}

char* test_determinant3(void){
    double a[9] = { 0.30112582,  0.67285109,  0.37710869,  
                    0.43101797,  0.73377317,  0.60827196,  
                    0.42625684,  0.3349345 ,  0.9300366 };
    double d = -0.014623566951869022;
    double t;

    t = determinant3( a );
    MU_ASSERT( "determinant OK", fabs( t - d) < 1e-6 );
    return 0;
}

char* test_inverse_mat3(void){
    double a[9] = { 0.30112582,  0.67285109,  0.37710869,  
                    0.43101797,  0.73377317,  0.60827196,  
                    0.42625684,  0.3349345 ,  0.9300366 };
    double e[9] = {-32.73514938,  34.15510178,  -9.06510735,
                  9.68179703,  -8.15894508,   1.41044748,
                 11.51655006, -12.71577263,   4.72202539 };
    double t[9];
    int i, ret;
    MU_ASSERT( "inverse runs", inverse_mat3( a, t) == 0 );
    ret = 0;
    for( i=0; i<9; i++){
        if (fabs(e[i]-t[i])/fabs(t[i])>1e-6){
            ret += 1;
        }
    }
    MU_ASSERT( "inverse correct", ret == 0 );
    return 0;
}


char* test_affine_transform_vec(void){
    double m[12] = { 0, 0, 1, 
                    0, 1, 0,
                    1, 0, 0,
                    0.1,0.2,0.3};
    double v[3] = { 3,4,5 };
    double r[3];
    affine_transform_vec( m, v, r );
    MU_ASSERT("transform_vec0",fabs(r[0] - 5.1)<1e-6);
    MU_ASSERT("transform_vec1",fabs(r[1] - 4.2)<1e-6);
    MU_ASSERT("transform_vec2",fabs(r[2] - 3.3)<1e-6);
    return 0;
}

char* test_mat3_prod(void){
    double a[9] = {  0.97300866,  0.53048923,  0.32333488,  
                    0.6900934 ,  0.44364110,  0.68283224,  
                    0.51383945,  0.87587622,  0.80127724 };
    double b[9] = {  0.06703767,  0.85926415,  0.38202935,  
                    0.68021704,  0.29180251,  0.38910406, 
                    0.59501276,  0.50158507,  0.07206758 };
    double c[9] = {  0.61846443,  1.15304949,  0.60143535, 
                    0.75432839,  1.06492656,  0.48546856,  
                    1.10700271,  1.09901539,  0.59485486 };
    double t[9];
    int i;
    mat3_prod( a, b, t);
    for( i=0 ; i<9 ; i++)
        MU_ASSERT( "mat3_prod", fabs((c[i]-t[i])/c[i])<1e-6 );
    return 0;
}


char* test_inverse_affine_mat(void){
    double a[3] = {4,5,6};
    double m[12] = { 0,0,1,  1,0,0, 0,1,0, 0.1, 0.2, 0.3 };
    double r[3];
    double im[12];
    double ir[3];
    int i;
    affine_transform_vec(m, a, r);
    inverse_affine_mat(m, im);
    affine_transform_vec(im, r, ir);
    for(i=0;i<3;i++)
        MU_ASSERT("inverse affine OK", fabs( ( ir[i] - a[i] )/a[i] ) < 1e-6);
    return 0;
}

char* test_norm3(void){
    double t[3] = {4,5,6};
    MU_ASSERT( "norm3 OK", fabs(norm3(t)- sqrt(77)) < 1e-6);
    return 0;
}

char* test_normalise_vector(void){
    double v[3] = {2,2,2};
    int i;
    normalise_vector(v);
    for(i=0;i<3;i++)
        MU_ASSERT("Normalise vector",fabs(3*v[i]*v[i] - 1)<1e-6);
    return 0;
}

char* test_mat3_axis_angle(void){
    double m[9];
    double a[3] = {0, 0, 1};
    double t[9],s30,c30;
    double p=30.0;
    int err,i;
    p = radians( 30.0 ) ;
    err = mat3_from_axis_angle( a, p, m );
    c30 = sqrt(3)/2.; /*  cos(30) */
    s30 = 0.5  ;      /*   sin(30) */
    t[0] = c30; t[1] = -s30; t[2]=0; 
    t[3] = s30; t[4] =  c30; t[5]=0; 
    t[6] = 0  ; t[7] =    0; t[8]=1;
    for(i=0;i<9;i++)
        MU_ASSERT("mat3 aa",fabs(m[i] - t[i])<1e-6);
    return 0;
}

char* test_rotate_vector_axis_angle(void){
    double a[9] = {0, 0, 1};
    MU_ASSERT("No axis angle test written",0);
    return 0;
}



char* test_suite(void){
    MU_RUN_TEST(test_mat3_transform_vec);
    MU_RUN_TEST(test_determinant3);
    MU_RUN_TEST(test_inverse_mat3);
    MU_RUN_TEST(test_affine_transform_vec);
    MU_RUN_TEST(test_mat3_prod);
    MU_RUN_TEST(test_inverse_affine_mat);
    MU_RUN_TEST(test_norm3);
    MU_RUN_TEST(test_normalise_vector);
    MU_RUN_TEST(test_mat3_axis_angle);
    MU_RUN_TEST(test_rotate_vector_axis_angle);
}

int main(int argc,
         char* argv[]){
    char *result = test_suite();
    printf("number of tests run: %d\n",tests_run);
    if( result ) printf("FAIL: %s\n", result);
    return 0!= result;
}

