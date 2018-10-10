#include <iostream>

using namespace std;

double calculate(double **arr, int rows, int cols) {
    int i, j, sum = 0, product;
    for(i = 0; i < rows; i++) {
        product = 1;
        for(j = 0; j < cols; j++)
            product *= arr[i][j];
      	  sum += product;
    }
    return sum;
}