#include <stdio.h>
#include <stdlib.h>

void FreeMatrix(int dim, double ***matrix) {
  for (int i = 0; i<dim; i++) {
    free((*matrix)[i]);
    (*matrix)[i] = NULL;
  }
  free(*matrix);
  *matrix = NULL;
}


//compute the inner product between two vectors
double InnerProduct(int dim, double *arr1, double *arr2) {
  double sum1 = 0;
  int i;
  for (i=0; i<dim; i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sum1;
}
