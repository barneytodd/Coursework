#include <stdio.h>
#include <stdlib.h>

//used to free the memory allocated for A and B when either an error arises or the end of the program is reached
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
