%module ddarray
%inline %{
// Helper function to create a 2d array


double **new_doubleddArrayTeste(int rows){
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

void delete_ddArray (double **arr, int rows, int cols){
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

double calculate(double **arr, int rows, int cols);

%}