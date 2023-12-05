#include <stdarg.h>
#include <stdio.h>
#include <math.h>

void GramSchmidt(int count, int dim, double A[][], ...) {
  va_list ap; //initialise list of variables
  int i, j, k; //initialise variables i, j, k
  va_start (ap, count); //initialise va_list
  for (i=0; i<count; i++) { //iterate through the variables (vectors)
    int *vector = va_arg (ap, int*); //store the vector in the variable vector
    for (k=0; k<dim; k++) { 
      A[i][k] = vector[k]; // initialise row i of A
      } 
    for (j=0; j<i; j++) { //iterate through the previous vectors
      double dot_product = 0.0; //initalise the dot product count
      for (k=0; k<dim; k++) {  //iterate through the entries in the vector to calculate the dot product
        dot_product += A[j][k] * vector[k]; //calculate the dot product
      } 
      for (k=0; k<dim; k++) {
        A[i][k] -= dot_product*A[j][k]; //subtract the dot_product times the jth normalised vector 
      }
    }
    //double norm = 0.0; //initialise the norm
    //for (j=0; j<dim; j++) {
    //  norm += A[i][j] * A[i][j]; //calculate the square sum of entries
    //}
    //norm = sqrt(norm); //sqrt to find norm
    //for (k=0; k<dim; k++) {
    //  A[i][k]/=norm; //normalise the entries in A
    //}
  } 
  va_end (ap);
}
  
    
void LLL(int count, double delta=3/4, int dim, double A[count][dim], ...) {
  va_list ap; //initialise list of variables
  int i, j, k; //initialise variables i, j, k
  va_start (ap, count); //initialise va_list
  for (i=0; i<count; i++) { //iterate through the variables (vectors)
    int *vector = va_arg (ap, int*); //store the vector in the variable vector
    for (k=0; k<dim; k++) { 
      A[i][k] = vector[k]; // initialise row i of A
      } 
  

int main() {
  double A[2][3];
  GramSchmidt(2, 3, A, (int[]){0, 0, 1}, (int[]){0, 2, 0});
  printf("Orthonormalized Vectors (A):\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
  return 0;
}
