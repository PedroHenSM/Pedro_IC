/* cec17.i */
%module cec17NoConstraints
%include "carrays.i" /* %array_functions(type,name) */
%array_functions (double, doubleArray);
%array_functions (long double, longDoubleArray);
/* Creates an void array | array_functions(type,name) */
/* %array_class(void,voidArray) Creates an void array for classes | array_class(type,name) */
/* %typemap(in) void* = double*; */

%{

/* Helper functions */

double **new_doubleddArray(int rows){
	double **arr = new double *[rows];
	return arr;
}

double **castToDouble(void *b){
	return (double**)b;
}

double **new_doubleddArray(int rows, int cols) {
    int i;
    double **arr = new double *[rows];
    for (i=0; i<rows; i++)
		arr[i] = new double[cols];
    return arr;
}

void delete_doubleddArray (double **arr, int rows, int cols){
	int i;
	for (i=0; i<rows; i++)
		delete[] arr[i];
	delete[] arr;
}

void doubleddArray_setitem(double **array, int row, int col, double value) {
    array[row][col] = value;
}

double doubleddArray_getitem(double **array, int row, int col) {
    return array[row][col];
}

void printDoubleArray(double *array, int size){
    int i;
    for (i = 0; i < size; i++){
        printf("%f ", array[i]);
    }
    printf("\n");
}

void printLongDoubleArray(long double *array, int size){
    int i;
    for (i = 0; i < size; i++){
        printf("%Lf ", array[i]);
    }
    printf("\n");
}

/* Includes headers files or function declarations */

// #include <WINDOWS.H>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

extern double *OShift,*M,*y,*z,*x_bound;
extern int ini_flag,n_flag,func_flag,*SS;

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

void sphere_func (double *, double *, int , double *,double *, int, int); /* Sphere */
void ellips_func(double *, double *, int , double *,double *, int, int); /* Ellipsoidal */
void bent_cigar_func(double *, double *, int , double *,double *, int, int); /* Discus */
void discus_func(double *, double *, int , double *,double *, int, int);  /* Bent_Cigar */
void dif_powers_func(double *, double *, int , double *,double *, int, int);  /* Different Powers */
void rosenbrock_func (double *, double *, int , double *,double *, int, int); /* Rosenbrock's */
void schaffer_F7_func (double *, double *, int , double *,double *, int, int); /* Schwefel's F7 */
void ackley_func (double *, double *, int , double *,double *, int, int); /* Ackley's */
void rastrigin_func (double *, double *, int , double *,double *, int, int); /* Rastrigin's  */
void weierstrass_func (double *, double *, int , double *,double *, int, int); /* Weierstrass's  */
void griewank_func (double *, double *, int , double *,double *, int, int); /* Griewank's  */
void schwefel_func (double *, double *, int , double *,double *, int, int); /* Schwefel's */
void katsuura_func (double *, double *, int , double *,double *, int, int); /* Katsuura */
void bi_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Lunacek Bi_rastrigin */
void grie_rosen_func (double *, double *, int , double *,double *, int, int); /* Griewank-Rosenbrock  */
void escaffer6_func (double *, double *, int , double *,double *, int, int); /* Expanded Scaffer¡¯s F6  */
void step_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Noncontinuous Rastrigin's  */
void happycat_func (double *, double *, int , double *,double *, int, int); /* HappyCat */
void hgbat_func (double *, double *, int , double *,double *, int, int); /* HGBat  */

/* New functions Noor Changes */
void sum_diff_pow_func(double *, double *, int , double *,double *, int, int); /* Sum of different power */
void zakharov_func(double *, double *, int , double *,double *, int, int); /* ZAKHAROV */
void levy_func(double *, double *, int , double *,double *, int, int); /* Levy */
void dixon_price_func(double *, double *, int , double *,double *, int, int); /* Dixon and Price */

void hf01 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 1 */
void hf02 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 2 */
void hf03 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 3 */
void hf04 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 4 */
void hf05 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 5 */
void hf06 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 6 */
void hf07 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 7 */
void hf08 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 8 */
void hf09 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 9 */
void hf10 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 10 */

void cf01 (double *, double *, int , double *,double *, int); /* Composition Function 1 */
void cf02 (double *, double *, int , double *,double *, int); /* Composition Function 2 */
void cf03 (double *, double *, int , double *,double *, int); /* Composition Function 3 */
void cf04 (double *, double *, int , double *,double *, int); /* Composition Function 4 */
void cf05 (double *, double *, int , double *,double *, int); /* Composition Function 5 */
void cf06 (double *, double *, int , double *,double *, int); /* Composition Function 6 */
void cf07 (double *, double *, int , double *,double *, int); /* Composition Function 7 */
void cf08 (double *, double *, int , double *,double *, int); /* Composition Function 8 */
void cf09 (double *, double *, int , double *,double *, int *, int); /* Composition Function 9 */
void cf10 (double *, double *, int , double *,double *, int *, int); /* Composition Function 10 */

void shiftfunc (double*,double*,int,double*);
void rotatefunc (double*,double*,int, double*);
void sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
void asyfunc (double *, double *x, int, double);
void oszfunc (double *, double *, int);
void cf_cal(double *, double *, int, double *,double *,double *,double *,int);

extern double *OShift,*M,*y,*z,*x_bound;
extern int ini_flag,n_flag,func_flag,*SS;

void cec17_test_func(double *, double *,int,int,int);


%}

/* Helper Functions */

double **new_doubleddArray(int rows);

double **castToDouble(void *b);

void delete_doubleddArray (double **arr, int rows, int cols);

void doubleddArray_setitem(double **array, int row, int col, double value);

double doubleddArray_getitem(double **array, int row, int col);

void printDoubleArray(double *array, int size);

void printLongDoubleArray(long double *array, int size);

// #include <WINDOWS.H>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

void sphere_func (double *, double *, int , double *,double *, int, int); /* Sphere */
void ellips_func(double *, double *, int , double *,double *, int, int); /* Ellipsoidal */
void bent_cigar_func(double *, double *, int , double *,double *, int, int); /* Discus */
void discus_func(double *, double *, int , double *,double *, int, int);  /* Bent_Cigar */
void dif_powers_func(double *, double *, int , double *,double *, int, int);  /* Different Powers */
void rosenbrock_func (double *, double *, int , double *,double *, int, int); /* Rosenbrock's */
void schaffer_F7_func (double *, double *, int , double *,double *, int, int); /* Schwefel's F7 */
void ackley_func (double *, double *, int , double *,double *, int, int); /* Ackley's */
void rastrigin_func (double *, double *, int , double *,double *, int, int); /* Rastrigin's  */
void weierstrass_func (double *, double *, int , double *,double *, int, int); /* Weierstrass's  */
void griewank_func (double *, double *, int , double *,double *, int, int); /* Griewank's  */
void schwefel_func (double *, double *, int , double *,double *, int, int); /* Schwefel's */
void katsuura_func (double *, double *, int , double *,double *, int, int); /* Katsuura */
void bi_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Lunacek Bi_rastrigin */
void grie_rosen_func (double *, double *, int , double *,double *, int, int); /* Griewank-Rosenbrock  */
void escaffer6_func (double *, double *, int , double *,double *, int, int); /* Expanded Scaffer¡¯s F6  */
void step_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Noncontinuous Rastrigin's  */
void happycat_func (double *, double *, int , double *,double *, int, int); /* HappyCat */
void hgbat_func (double *, double *, int , double *,double *, int, int); /* HGBat  */

/* New functions Noor Changes */
void sum_diff_pow_func(double *, double *, int , double *,double *, int, int); /* Sum of different power */
void zakharov_func(double *, double *, int , double *,double *, int, int); /* ZAKHAROV */
void levy_func(double *, double *, int , double *,double *, int, int); /* Levy */
void dixon_price_func(double *, double *, int , double *,double *, int, int); /* Dixon and Price */

void hf01 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 1 */
void hf02 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 2 */
void hf03 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 3 */
void hf04 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 4 */
void hf05 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 5 */
void hf06 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 6 */
void hf07 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 7 */
void hf08 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 8 */
void hf09 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 9 */
void hf10 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 10 */

void cf01 (double *, double *, int , double *,double *, int); /* Composition Function 1 */
void cf02 (double *, double *, int , double *,double *, int); /* Composition Function 2 */
void cf03 (double *, double *, int , double *,double *, int); /* Composition Function 3 */
void cf04 (double *, double *, int , double *,double *, int); /* Composition Function 4 */
void cf05 (double *, double *, int , double *,double *, int); /* Composition Function 5 */
void cf06 (double *, double *, int , double *,double *, int); /* Composition Function 6 */
void cf07 (double *, double *, int , double *,double *, int); /* Composition Function 7 */
void cf08 (double *, double *, int , double *,double *, int); /* Composition Function 8 */
void cf09 (double *, double *, int , double *,double *, int *, int); /* Composition Function 9 */
void cf10 (double *, double *, int , double *,double *, int *, int); /* Composition Function 10 */

void shiftfunc (double*,double*,int,double*);
void rotatefunc (double*,double*,int, double*);
void sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
void asyfunc (double *, double *x, int, double);
void oszfunc (double *, double *, int);
void cf_cal(double *, double *, int, double *,double *,double *,double *,int);

void cec17_test_func(double *, double *,int,int,int);

