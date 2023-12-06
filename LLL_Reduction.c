#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "LLL_Reduction.h"

void GramSchmidt(int count, int dim, double A[][count]) {
  int i, j, k; //initialise variables i, j, k
  double B[dim][count]
  for (i=0; i<count; i++) { //iterate through the variables (vectors)
    for (j=0; j<i; j++) { //iterate through the previous vectors
      double inner_product = 0.0; //initalise the dot product count
      for (k=0; k<dim; k++) {  //iterate through the entries in the vector to calculate the dot product
        inner_product += A[k][j] * A[k][i]; //calculate the dot product
      } 
      for (k=0; k<dim; k++) {
        A[k][i] -= inner_product*A[k][j]; //subtract the dot_product times the jth normalised vector 
      }
    }
    //double norm = 0.0; //initialise the norm
    //for (j=0; j<dim; j++) {
    //  norm += A[j][i] * A[j][i]; //calculate the square sum of entries
    //}
    //norm = sqrt(norm); //sqrt to find norm
    //for (k=0; k<dim; k++) {
    //  A[k][i]/=norm; //normalise the entries in A
    //}
  } 
}
  

void update_matrices(int count, int dim, double A[][count], double B[][count], double M[][count]) {
  int i, j, k;
  for (i=0; i<count; i++) {
    for (j=0; j<dim; j++) {
      B[j][i] = A[j][i];  
    }
  }
  GramSchmidt(count, dim, B);
  for (i=0; i<count; i++) { 
    for (j=0; j<count; j++) {
      double inner_product1 = 0.0;
      double inner_product2 = 0.0;
      for (k=0; k<dim; k++) {
        inner_product1 += A[k][i] * B[k][j];
        inner_product2 += B[k][j] * B[k][j];
      }
      M[i][j] = inner_product1/inner_product2; // calculate components of M
    }
  }
}

void LLL(int count, double delta, int dim, double A[][dim], double B[][dim], ...) {
  va_list ap; //initialise list of variables
  int i, j, k, m; //initialise variables i, j, k
  va_start (ap, B); //initialise va_list
  for (i=0; i<count; i++) { //iterate through the variables (vectors)
    printf("i: %d\n", i);
    int *vector = va_arg (ap, int*); //store the vector in the variable vector
    for (k=0; k<dim; k++) { 
      printf("k: %d\n", k);
      A[k][i] = vector[k]; // initialise row i of A
      B[k][i] = vector[k]; // initialise to be the same as A
      printf("Vector Entries: %d\t", vector[k]);
      } 
  }
  printf("A, B after initialisation:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("A: %.4f\t", A[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("B: %.4f\t", B[i][j]);
        }
        printf("\n");
    }
  GramSchmidt(count, dim, B); //GramSchmidt B
  printf("A, B after GramSchmidt:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("A: %.4f\t", A[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("B: %.4f\t", B[i][j]);
        }
        printf("\n");
    }
  double M[count][count]; // initialise a new matrix M
  for (i=0; i<count; i++) { 
    for (j=0; j<count; j++) {
      double inner_product1 = 0.0;
      double inner_product2 = 0.0;
      for (k=0; k<dim; k++) {
        inner_product1 += A[k][i] * B[k][j];
        inner_product2 += B[k][j] * B[k][j];
      }
      M[i][j] = inner_product1/inner_product2; // calculate components of M
    }
  }
  printf("M after initialisation:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.4f\t", M[i][j]);
        }
        printf("\n");
    }
  k = 2;
  while (k<=dim) {
    for (j=k-1; j>0; j--) {
      if (fabs(M[k][j]) > 1/2) {
        for (i=0; i<dim; i++) {
          A[i][k] -= M[k][j] * A[i][j];
        }
        update_matrices(count, dim, A, B, M);
      }
    }
    printf("A, B, M after updating:\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("A: %.4f\t", A[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("B: %.4f\t", B[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("M: %.4f\t", M[i][j]);
        }
        printf("\n");
    }
    double inner_product1 = 0.0;
    double inner_product2 = 0.0;
    for (i=0; i<dim; i++) {
      inner_product1 += B[i][k]*B[i][k];
      inner_product2 += B[i][k-1]*B[i][k-1];
      }
    printf("Inner product 1: %.4f\n", inner_product1);
    printf("Compare to: %.4f\n", ((delta - (M[k][k-1]*M[k][k-1])) * inner_product2));
    if (inner_product1 > ((delta - (M[k][k-1]*M[k][k-1])) * inner_product2)) {
      k+=1;
        }
      else {
        double temp[dim];
        for (i=0; i<dim; i++) {
          temp[i] = A[i][k];
          A[i][k] = A[i][k-1];
          A[i][k-1] = temp[i];
        }
        update_matrices(count, dim, A, B, M); 
        printf("A, B, M after updating again:\n");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                printf("A: %.4f\t", A[i][j]);
            }
            printf("\n");
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                printf("B: %.4f\t", B[i][j]);
            }
            printf("\n");
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                printf("M: %.4f\t", M[i][j]);
            }
            printf("\n");
        }
      }
  }
  va_end(ap);
}
  
  
  


