#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "LLL_Reduction.h"

//compute the inner product between two vectors
double InnerProduct(int dim, double *arr1, double *arr2) {
  double sum1 = 0;
  int i;
  for (i=0; i<dim; i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sum1;
}

//compute GramSchmidt orthogonalisation without normalisation
void GramSchmidt(int dim, int start, double B[][dim]) {
  int i, j, k; 
  double mu_ij;
  double vec1[dim]; //store values to subtract from initial vectors
  
  //iterate through the initial vectors
  for (i=start; i<dim; i++) { 
    for (j=0; j<dim; j++) {
      vec1[j] = 0;
    }
    //iterate through the previous vectors
    for (j=0; j<i; j++) { 
      mu_ij = InnerProduct(dim, B[i], B[j])/InnerProduct(dim, B[j], B[j]);
      for (k=0; k<dim; k++) {
        vec1[k] += mu_ij * B[j][k]; //add the dot_product times the jth normalised vector 
      }
      
    }
    //subtract from the ith initial vector
    for (k=0; k<dim; k++) {
      B[i][k] -= vec1[k];
    }
  }
}
  
//when A gets updated, recompute B to be the GramSchmidt orthogonalised version of the updated A
void update_matrices(int dim, int start, double **A, double B[][dim]) {
  int i, j;

  //set B to equal A
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      B[i][j] = A[i][j];  
    }
  }
  
  GramSchmidt(dim, start, B);
}


/// Lenstra–Lenstra–Lovász reduce the input matrix A
void LLL(double delta, int dim, double **A) {
  
  int i, j, k; //initialise variables i, j, k
  double B[dim][dim];
  
  //set be to be equal to A
  for (i=0; i<dim; i++) { 
    for (j=0; j<dim; j++) { 
      B[i][j] = A[i][j]; 
      } 
  }

  //GramSchmidt orthogonslise B
  GramSchmidt(dim, 0, B); 
  
  k = 1;
  int m = 0;
  double mu_kj;
  double mu_k_kminus1; 

  //
  while (k<dim) {
    //reduce the kth vector until for all j<k, mu_kj<=0.5
    for (j=k-1; j>=0; j--) {
      mu_kj = InnerProduct(dim, A[k], B[j])/InnerProduct(dim, B[j], B[j]); 
      if (fabs(mu_kj) > 0.5) {
        for (i=0; i<dim; i++) {
          A[k][i] -= round(mu_kj) * A[j][i];
        }
        
        update_matrices(dim, k, A, B);
      }
    }

    //if 
    mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]); 
    if (InnerProduct(dim, B[k], B[k]) > ((delta - (mu_k_kminus1*mu_k_kminus1)) * InnerProduct(dim, B[k-1], B[k-1]))) {
      k+=1;
        }
    else {
      //swap A[k] and A[k-1]
      for (i=0; i<dim; i++) {
        A[k][i] += A[k-1][i];
        A[k-1][i] = A[k][i] - A[k-1][i];
        A[k][i] -= A[k-1][i];
      }
      update_matrices(dim, k, A, B); 
      k = fmax(k-1, 1);
      //mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]);           
    }
    m+=1;
    if (m==10) {
      break;
    }      
  }
}
  
  
  


