#include <stdarg.h> //check which ones we need
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "LLL_Reduction.h"

double InnerProduct(int dim, double *arr1, double *arr2) {
  double sum1 = 0;
  int i;
  for (i=0; i<dim; i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sum1;
}

void GramSchmidt(int dim, int start, double B[][dim]) {
  int i, j, k; //initialise variables i, j, k
  double mu_ij;
  double vec1[dim];
  for (i=0; i<dim; i++) {
    vec1[i] = 0;
  }
  for (i=start; i<dim; i++) { //iterate through the variables (vectors)
    for (j=0; j<i; j++) { //iterate through the previous vectors
      mu_ij = InnerProduct(dim, B[i], B[j])/InnerProduct(dim, B[j], B[j]);
      for (k=0; k<dim; k++) {
      }
      for (k=0; k<dim; k++) {
        vec1[k] += mu_ij * B[j][k]; //subtract the dot_product times the jth normalised vector 
      }
    }
    for (k=0; k<dim; k++) {
      B[i][k] -= vec1[k];
    }
  }
}
  

void update_matrices(int dim, int start, double (*A)[dim], double B[][dim]) {
  int i, j;
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      B[i][j] = A[i][j];  
    }
  }
  GramSchmidt(dim, start, B);
}


//change to produce row vectors
void LLL(double delta, int dim, double (*A)[dim], ...) {
  va_list ap; //initialise list of variables
  int i, j, k; //initialise variables i, j, k
  double B[dim][dim];
  double *vector;
  va_start (ap, A); //initialise va_list
  for (i=0; i<dim; i++) { //iterate through the variables (vectors)
    vector = va_arg (ap, double *); //store the vector in the variable vector
    for (k=0; k<dim; k++) { 
      A[i][k] = vector[k]; // initialise row i of A
      B[i][k] = vector[k]; // initialise to be the same as A
      } 
  }
  printf("A after initialisation:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("A: %.4f\t", A[i][j]);
        }
        printf("\n");
    }
  GramSchmidt(dim, 0, B); //GramSchmidt B
  
  k = 1;
  int m = 0;
  double mu_kj;
  double mu_k_kminus1; //= InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]);
  while (k<dim) {
    printf("k: %d\n", k);
    for (j=k-1; j>=0; j--) {
      mu_kj = InnerProduct(dim, A[k], B[j])/InnerProduct(dim, B[j], B[j]); 
      if (fabs(mu_kj) > 0.5) {
        for (i=0; i<dim; i++) {
          printf("A[k][i] before: %.4f\n", A[k][i]);
          A[k][i] -= round(mu_kj) * A[j][i];
          printf("A[k][i] after: %.4f\n", A[k][i]);
        }
        printf("A after updating:\n");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                printf("A: %.4f\t", A[i][j]);
            }
            printf("\n");
        }
        update_matrices(dim, k, A, B);
      }
    }
    
    
    
    mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]); 
    printf("compare1: %.4f\n", InnerProduct(dim, B[k], B[k]));
    printf("compare2: %.4f\n", ((delta - (mu_k_kminus1*mu_k_kminus1)) * InnerProduct(dim, B[k-1], B[k-1])));
    if (InnerProduct(dim, B[k], B[k]) > ((delta - (mu_k_kminus1*mu_k_kminus1)) * InnerProduct(dim, B[k-1], B[k-1]))) {
      printf("yes\n");
      k+=1;
      printf("k2: %d", k);
      mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]); 
        }
    else {
        A[k][i] += A[k-1][i];
        A[k-1][i] = A[k][i] - A[k-1][i];
        A[k][i] -= A[k-1][i];
      }
      update_matrices(dim, k, A, B); 
      printf("A after updating again:\n");
      for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
              printf("A: %.4f\t", A[i][j]);
          }
          printf("\n");
      k = fmax(k-1, 1);
      mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]);           
    }
    m+=1;
    if (m==10) {
      break;
    }      
  }
  va_end(ap);
}
  
  
  


