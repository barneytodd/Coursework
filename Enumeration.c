#include <string.h>

double VectorNorm(int dim, int *arr) {
  double sum1;
  int i;
  for (i=0, i<dim, i++) {
    sum1 += arr[i]*arr[i];
  }
  return sqrt(sum1);
}

double ShortestVector(int dim, double (*A)[dim]) { //A is  amtrix of row vectors stacked
  int i, j, k, l;
  bool boolarray_prev[dim]; //checks if each vector was changed in the previous loop
  bool boolarray_curr[dim]; //checks if each vector has been changed in the current loop
  for (i=0; i<dim; i++) {
    boolarray_prev[i] = True; //initialises boolarray to all True
    boolarray_curr[i] = False;
  }
  double temp_vector[dim];
  double i_length
  double j_length
  double temp_length
  for (i=0; i<dim; i++) {
    for (j=i; j<dim; j++) {
      if (i!=j && boolarray_prev[j]) {
        i_length = VectorNorm(dim, A[i]);
        j_lenght = VectorNorm(dim, A[j]);
        for (k=0; k<dim; k++) {
          temp_vector[k] = A[i][k] + A[j][k];
          temp_length = VectorNorm(dim, temp_vector); 
          if (temp_length < i_length) {
            memcpy(A[i], temp_vector, dim*sizeof(double));
            boolarray_prev[i] = boolarray_curr[i] = True;
          }
          else if (temp_length < j_length) {
            memcpy(A[j], temp_vector, dim*sizeof(double));
            boolarray_prev[j] = boolarray_curr[j] = True;
          }
          else {
            temp_vector[k] = A[i][k] - A[j][k];
            temp_length = VectorNorm(dim, temp_vector); 
            if (temp_length < i_length) {
            memcpy(A[i], temp_vector, dim*sizeof(double));
            boolarray_prev[i] = boolarray_curr[i] = True;
          }
          else if (temp_length < j_length) {
            memcpy(A[j], temp_vector, dim*sizeof(double));
            boolarray_prev[j] = boolarray_curr[j] = True;
          }
          }
        }
      }
    }
  }
  memcpy(boolarray_prev, boolarray_curr, dim*sizeof(bool));
}
  
