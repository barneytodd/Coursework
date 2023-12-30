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
  
  for (i=start; i<dim; i++) { //iterate through the variables (vectors)
    for (j=0; j<dim; j++) {
      vec1[j] = 0;
    }
    for (j=0; j<i; j++) { //iterate through the previous vectors
      mu_ij = InnerProduct(dim, B[i], B[j])/InnerProduct(dim, B[j], B[j]);
      //
      //if (j==0) {
      //  printf("mu_ij: %.4f %d %d\n", mu_ij, i, j);
      //}
      //for (k=0; k<dim; k++) {
      //}
      for (k=0; k<dim; k++) {
        if (k==0 && mu_ij == (B[i][0] * B[j][0])/(B[j][0] * B[j][0])) {
          //vec[k] += 0;
          continue;
        }
        vec1[k] += mu_ij * B[j][k]; //subtract the dot_product times the jth normalised vector 
      }
      printf("vec1\n");
      for (k=0; k<dim; k++) {
        printf("%.4f\t", vec1[k]);
      }
      printf("\n");
          
          //if (j==0 || j==1) {
            //printf("xyz %.4f\n", (B[i][0] * B[j][0])/(B[j][0] * B[j][0]));
            //printf("xyz %.4f\n",  mu_ij * B[j][0] - B[i][0]);
            //printf("xyz %.4f\n",  * B[j][k] - B[i][0]);
            //printf("j, vec1[0]+ %d, %.4f\n", j, mu_ij * B[j][k]);
          //}
        
      
      
    }
    //printf("vec1[0](2) %.4f\n", vec1[0]);
    //printf("%.4f\n", B[i][0]);
    //for (k=0; k<dim; k++) {
    //  B[i][k] -= vec1[k];
    //}
    //printf("B[i]\n");
    //for (k=0;k<dim;k++) {
    //  printf("%.4f %d %d\t", B[i][k], i, k);
    //}
    //printf("\n");
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
void LLL(double delta, int dim, double (*A)[dim]) {
  //va_list ap; //initialise list of variables
  int i, j, k; //initialise variables i, j, k
  double B[dim][dim];
  //double *vector;
  //va_start (ap, A); //initialise va_list
  for (i=0; i<dim; i++) { //iterate through the variables (vectors)
    //vector = va_arg (ap, double *); //store the vector in the variable vector
    for (j=0; j<dim; j++) { 
      //A[i][k] = vector[k]; // initialise row i of A
      B[i][j] = A[i][j]; //vector[k]; // initialise to be the same as A
      } 
  }
  GramSchmidt(dim, 0, B); //GramSchmidt B

  //printf("B1:\n");
  //  for (int i = 0; i < dim; i++) {
  //      for (int j = 0; j < dim; j++) {
  //          printf("%.4f\t", B[i][j]);
  //      }
  //      printf("\n");
  //  }


  
  k = 1;
  int m = 0;
  double mu_kj;
  double mu_k_kminus1; //= InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]);
  while (k<dim) {
    for (j=k-1; j>=0; j--) {
      mu_kj = InnerProduct(dim, A[k], B[j])/InnerProduct(dim, B[j], B[j]); 
      if (fabs(mu_kj) > 0.5) {
        for (i=0; i<dim; i++) {
          A[k][i] -= round(mu_kj) * A[j][i];
        }
        update_matrices(dim, k, A, B);
      }
    }

    //printf("After while (A):\n");
    //for (int i = 0; i < dim; i++) {
    //    for (int j = 0; j < dim; j++) {
    //        printf("%.4f\t", A[i][j]);
    //    }
    //    printf("\n");
    //}

    //printf("After while (B):\n");
    //for (int i = 0; i < dim; i++) {
    //    for (int j = 0; j < dim; j++) {
    //        printf("%.4f\t", A[i][j]);
    //    }
    //    printf("\n");
    //}
    
    
    mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]); 
    if (InnerProduct(dim, B[k], B[k]) > ((delta - (mu_k_kminus1*mu_k_kminus1)) * InnerProduct(dim, B[k-1], B[k-1]))) {
      k+=1;
      mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]); 
        }
    else {
      for (i=0; i<dim; i++) {
        A[k][i] += A[k-1][i];
        A[k-1][i] = A[k][i] - A[k-1][i];
        A[k][i] -= A[k-1][i];
      }
      update_matrices(dim, k, A, B); 
      k = fmax(k-1, 1);
      mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]);           
    }
    m+=1;
    if (m==10) {
      break;
    }      
  }
  //va_end(ap);
}
  
  
  


