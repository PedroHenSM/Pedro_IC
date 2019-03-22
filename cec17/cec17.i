/* cec17.i */
%module cec17
%include "carrays.i" /* %array_functions(type,name) */
%array_functions (double, doubleArray);
/* Creates an void array | array_functions(type,name) */
/* %array_class(void,voidArray) Creates an void array for classes | array_class(type,name) */
/* %typemap(in) void* = double*; */
%{
/* Includes headers files or function declarations */

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <iostream>
using namespace std;

extern double *OShift,*M,*M1,*M2,*y,*z,*z1,*z2;
extern int ini_flag,n_flag,func_flag,f5_flag;

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029



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

/* CEC17 */

void cec17_test_COP(double *x, double *f, double *g,double *h, int nx, int mx,int func_num);
void COP_01 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O1 */
void COP_02 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O2 */
void COP_03 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O3 */
void COP_04 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O4 */
void COP_05 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O5 */
void COP_06 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O6 */
void COP_07 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O7 */
void COP_08 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O8 */
void COP_09 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O9 */
void COP_10 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_10 */
void COP_11 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_11 */
void COP_12 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_12 */
void COP_13 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_13 */
void COP_14 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_14 */
void COP_15 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_15 */
void COP_16 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_15 */
void COP_17 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_17 */
void COP_18 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_18 */
void COP_19 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_19 */
void COP_20 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_20 */
/* Note that COP21-COP28 are  the rotated versioins of COP12-COP18. */
void loadShiftData(int fun_num, int dim, double *pV );
void loadRotateData(int func_num, int dim, double *pM );
void shiftfunc (double*,double*,int,double*);
void rotatefunc (double*,double*,int, double*);
void sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
int  sgn(double val);
double round(double val);


%}

// Helper function to create a 2d array

double **new_doubleddArray(int rows);

double **castToDouble(void *b);

void delete_doubleddArray (double **arr, int rows, int cols);

void doubleddArray_setitem(double **array, int row, int col, double value);

double doubleddArray_getitem(double **array, int row, int col);

void cec17_test_COP(double *x, double *f, double *g,double *h, int nx, int mx,int func_num);
void COP_01 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O1 */
void COP_02 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O2 */
void COP_03 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O3 */
void COP_04 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O4 */
void COP_05 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O5 */
void COP_06 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O6 */
void COP_07 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O7 */
void COP_08 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O8 */
void COP_09 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_O9 */
void COP_10 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_10 */
void COP_11 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_11 */
void COP_12 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_12 */
void COP_13 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_13 */
void COP_14 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_14 */
void COP_15 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_15 */
void COP_16 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_15 */
void COP_17 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_17 */
void COP_18 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_18 */
void COP_19 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_19 */
void COP_20 (double *, double *, double *, double *, int , double *,double *, int, int); /* COP_20 */
/* Note that COP21-COP28 are  the rotated versioins of COP12-COP18. */
void loadShiftData(int fun_num, int dim, double *pV );
void loadRotateData(int func_num, int dim, double *pM );
void shiftfunc (double*,double*,int,double*);
void rotatefunc (double*,double*,int, double*);
void sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
int  sgn(double val);
double round(double val);

