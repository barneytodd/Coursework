#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "LLL_Reduction.h"
#include <stdbool.h>

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
  double mag1;
  double mag2;
  //iterate through the initial vectors
  for (i=start; i<dim; i++) { 
    mag1 = sqrt(InnerProduct(dim, B[i], B[i]));
    for (j=0; j<dim; j++) {
      vec1[j] = 0;
      B[i][j] /= mag1;
    }
    //iterate through the previous vectors
    for (j=0; j<i; j++) { 
      mag2 = sqrt(InnerProduct(dim, B[j], B[j]));
      for (k=0;k<dim;k++) {
        B[j][k] /= mag2;
      }
      if (j<2 && i<5) {
        printf("j: %d\t", j);
        printf("ip check3: %.40f\n", InnerProduct(dim, B[i], B[j]) / InnerProduct(dim, B[j], B[j]) - B[i][0]/B[j][0]);
      }
      mu_ij = InnerProduct(dim, B[i], B[j])/InnerProduct(dim, B[j], B[j]);
      for (k=0; k<dim; k++) {
        vec1[k] += mu_ij * B[j][k] * mag1 / mag2; //add the dot_product times the jth normalised vector 
      }
      if (j==0 && i<5) {
        printf("i: %d\t", i);
        printf("vec1[0]: %.4f\t", vec1[0]);
        printf("B[i][0]: %.4f\n", B[i][0]);
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
  //for (i=0; i<dim; i++) { 
  //  for (j=0; j<dim; j++) { 
  //    B[i][j] = A[i][j]; 
  //    } 
  //}

  //GramSchmidt orthogonslise B
  //GramSchmidt(dim, 0, B); 
  update_matrices(dim, 0, A, B);
  printf("B: \n");
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      printf("%.4f\t", B[i][j]);
    }
    printf("\n");
  }
  
  k = 1;
  double mu_kj;
  double mu_k_kminus1; 
  bool zero_check = false; //checks if any vectors are reduced to 0, this can happen if the input vectors are linearly dependent and can result in an infinite loop
  //iterate through the LLL Reduction steps until:
  //(B[k] . B[k]) > (delta - mu_k_k-1) * (B[k-1] . B[k-1]) for every k, and
  //mu_kj<=0.5 for all k, j<k
  int m = 0;
  while (k<dim) {
    //reduce the kth vector until for all j<k, mu_kj<=0.5
    for (j=k-1; j>=0; j--) {
      mu_kj = InnerProduct(dim, A[k], B[j])/InnerProduct(dim, B[j], B[j]); 
      if (fabs(mu_kj) > 0.5) {
        zero_check = true;
        for (i=0; i<dim; i++) {
          A[k][i] -= round(mu_kj) * A[j][i];
          if (A[k][i] != 0) zero_check = false;
        }
        if (zero_check == true) {
          printf("Error: input vectors are linearly dependent\n");
          exit(1);
        }
        update_matrices(dim, k, A, B);
      }
    }
    //LLL basis reduction requires (B[k] . B[k]) > (delta - mu_k_k-1) * (B[k-1] . B[k-1]) for every k
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
      update_matrices(dim, k-1, A, B); 
      k = fmax(k-1, 1);          
    }
    m++;
    if (m % 10 == 0) { //need to improve this
      printf("While loop failed\n");
        exit(1);
    }
  }
  printf("m: %d\n", m);
}
  
  
  


