#include <stdarg.h> //check which ones we need
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "LLL_Reduction.h"

double InnerProduct(int dim, double *arr1, double *arr2) {
  //int k;
  //for (k=0; k<dim; k++) {
  //      printf("B[i]: %.4f\t", arr1[k]);
  //      printf("\n");
  //      printf("B[j]: %.4f\t", arr2[k]);
  //      printf("\n");
  //}
  double sum1 = 0;
  int i;
  for (i=0; i<dim; i++) {
    //printf("arr1[i]: %.4f\n", arr1[i]);
    //printf("arr2[i]: %.4f\n", arr2[i]);
    //printf("product: %.4f\n", arr1[i]*arr2[i]);
    sum1 += arr1[i]*arr2[i];
  }
  printf("sum1: %.4f\n", sum1);
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
      printf("mu_ij: %.4f\n", mu_ij);
      //printf("ip1: %.4f\n", InnerProduct(dim, B[i], B[j]));
      //printf("ip2: %.4f\n", InnerProduct(dim, B[j], B[j]));
      for (k=0; k<dim; k++) {
        //printf("B[i]: %.4f\t", B[i][k]);
        //printf("\n");
        //printf("B[j]: %.4f\t", B[j][k]);
        //printf("\n");
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
  //printf("updating before GS:\n");
    //for (int i = 0; i < 3; i++) {
      //  for (int j = 0; j < 3; j++) {
        //    printf("B: %.4f\t", B[i][j]);
        //}
        //printf("\n");
    //}
  GramSchmidt(dim, start, B);
  //printf("updating after GS:\n");
    //for (int i = 0; i < 3; i++) {
      //  for (int j = 0; j < 3; j++) {
        //    printf("B: %.4f\t", B[i][j]);
        //}
        //printf("\n");
    //}
  //for (i=0; i<count; i++) { 
  //  for (j=0; j<count; j++) {
  //    double inner_product1 = 0.0;
  //    double inner_product2 = 0.0;
  //    for (k=0; k<dim; k++) {
  //      inner_product1 += A[k][i] * B[k][j];
  //      inner_product2 += B[k][j] * B[k][j];
  //    }
  //    M[i][j] = inner_product1/inner_product2; // calculate components of M
  //  }
  //}
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
  GramSchmidt(dim, 0, B); //GramSchmidt B
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
  //double M[count][count]; // initialise a new matrix M
  //for (i=0; i<count; i++) { 
    //for (j=0; j<count; j++) {
      //double inner_product1 = 0.0;
      //double inner_product2 = 0.0;
      //for (k=0; k<dim; k++) {
        //inner_product1 += A[k][i] * B[k][j];
        //inner_product2 += B[k][j] * B[k][j];
      //}
      //M[i][j] = inner_product1/inner_product2; // calculate components of M
    //}
  //}
  //printf("M after initialisation:\n");
    //for (int i = 0; i < 3; i++) {
        //for (int j = 0; j < 3; j++) {
            //printf("%.4f\t", M[i][j]);
        //}
        //printf("\n");
    //}
  k = 1;
  int m = 0;
  double mu_kj;
  double mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]);
  while (k<dim) {
    for (j=k-1; j>=0; j--) {
      mu_kj = InnerProduct(dim, A[k], B[j])/InnerProduct(dim, B[j], B[j]); 
      printf("k: %d\n", k);
      printf("j: %d\n", j);
      printf("A, B:\n");
      for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
              printf("A: %.4f\t", A[i][j]);
          }
          printf("\n");
      }
      //printf("A, B after initialisation:\n");
      for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
              printf("B: %.4f\t", B[i][j]);
          }
          printf("\n");
      }
      printf("mu_kj: %.4f\n", mu_kj);
      if (fabs(mu_kj) > 1/2) {
        //int Mint = round(M[k][j]);
        for (i=0; i<dim; i++) {
          A[k][i] -= round(mu_kj) * A[j][i];
        }
        update_matrices(dim, k, A, B);
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
    //for (int i = 0; i < 3; i++) {
      //  for (int j = 0; j < 3; j++) {
        //    printf("M: %.4f\t", M[i][j]);
        //}
        //printf("\n");
    //}
    //double inner_product1 = 0.0;
    //double inner_product2 = 0.0;
    //for (i=0; i<dim; i++) {
      //inner_product1 += B[i][k]*B[i][k];
      //printf("B[i][k]: %.4f, i: %d, k: %d\n", B[i][k], i, k);
      //inner_product2 += B[i][k-1]*B[i][k-1];
      //}
    //printf("Inner product 1: %.4f\n", inner_product1);
    //printf("Compare to: %.4f\n", ((delta - (M[k][k-1]*M[k][k-1])) * inner_product2));
    //mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]); 
    printf("compare1: %.4f\n", InnerProduct(dim, B[k], B[k]));
    printf("compare2: %.4f\n", ((delta - (mu_k_kminus1*mu_k_kminus1)) * InnerProduct(dim, B[k-1], B[k-1])));
    printf("k: %d\n", k);
    if (InnerProduct(dim, B[k], B[k]) > ((delta - (mu_k_kminus1*mu_k_kminus1)) * InnerProduct(dim, B[k-1], B[k-1]))) {
      k+=1;
      mu_k_kminus1 = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]); 
        }
    else {
      for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
          printf("A[k]: %.4f\t", A[k][j]);
          printf("\n");
        }
        for (j=0; j<dim; j++) {
          printf("A[k-1]: %.4f\t", A[k-1][j]);
          printf("\n");
        }
        A[k][i] += A[k-1][i];
        A[k-1][i] = A[k][i] - A[k-1][i];
        A[k][i] -= A[k-1][i];
      }
      update_matrices(dim, k, A, B); 
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
      //for (int i = 0; i < 3; i++) {
        //  for (int j = 0; j < 3; j++) {
          //    printf("M: %.4f\t", M[i][j]);
          //}
          //printf("\n");
      //}
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
  
  
  


