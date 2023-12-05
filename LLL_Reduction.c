#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void GramSchmidt(int count, int dim, double A[][]) {
  int i, j, k; //initialise variables i, j, k
  for (i=0; i<count; i++) { //iterate through the variables (vectors)
    for (j=0; j<i; j++) { //iterate through the previous vectors
      double inner_product = 0.0; //initalise the dot product count
      for (k=0; k<dim; k++) {  //iterate through the entries in the vector to calculate the dot product
        inner_product += A[j][k] * A[i][k]; //calculate the dot product
      } 
      for (k=0; k<dim; k++) {
        A[i][k] -= inner_product*A[j][k]; //subtract the dot_product times the jth normalised vector 
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
  
    
void LLL(int count, double delta=3/4, int dim, double A[count][dim], double B[count][dim], ...) {
  va_list ap; //initialise list of variables
  int i, j, k, m; //initialise variables i, j, k
  va_start (ap, count); //initialise va_list
  for (i=0; i<count; i++) { //iterate through the variables (vectors)
    int *vector = va_arg (ap, int*); //store the vector in the variable vector
    for (k=0; k<dim; k++) { 
      A[i][k] = vector[k]; // initialise row i of A
      B[i][k] = vector[k]; // initialise to be the same as A
      } 
  }
  GramShmidt(count, dim, B); //GramSchmidt B
  double M[count][dim]; // initialise a new matrix M
  for (i=0; i<count; i++) { 
    for (j=0; j<count; j++) {
      double inner_product1 = 0.0;
      double inner_product2 = 0.0;
      for (k=0; k<dim; k++) {
        inner_product1 += A[i][k] * B[j][k];
        inner_product2 += B[j][k] * B[j][k];
      }
      M[i][j] = inner_product1/inner_product2; // calculate components of M
    }
  }
  k = 2;
  while (k<=dim) {
    for (j=k-1, j>0, j--) {
      if (abs(M[k][j]) > 1/2) {
        for (i=0; i<dim; i++) {
          A[k][i] -= M[k][j] * A[j][i]
        }
        for (i=0,; i<dim; i++) {
          for (m=0; m<count; m++) {
            B[m][i] = A[m][i];  
          }    
        }
        GramSchmidt(count, dim, B);
        //need to update M
      }
    }
    double inner_product1 = 0.0;
    double inner_product2 = 0.0;
    for (i=0; i<dim; i++) {
      inner_product1 += B[k][i]*B[k][i];
      inner_product2 += B[k-1][i]*B[k-1][i];
      }
    if (inner_product1 > ((delta - (M[k][k-1]*M[k][k-1])) * inner_product2) {
      k+=1
    }
  }
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
