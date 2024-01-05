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

bool CheckOrth(int dim, int start, double **B) {
  int i, j, k;
  double mag1;
  double mag2;
  for (i=fmax(start, 1); i<dim; i++) {
      mag1 = sqrt(InnerProduct(dim, B[i], B[i]));
      for (k=0; k<dim; k++) {
        B[i][k] /= mag1;
      }
      for (j=0; j<i; j++) {
        mag2 = sqrt(InnerProduct(dim, B[j], B[j]));
        for (k=0; k<dim; k++) {
          B[j][k] /= mag2;
        }
        if (fabs(InnerProduct(dim, B[i], B[j])) > 0.01) {
          printf("Failed: %.4f\n", fabs(InnerProduct(dim, B[i], B[j])));
					for (k=0; k<dim; k++) {
						B[j][k] *= mag2;
						B[i][k] *= mag1;
					}
					return false;
        }
        for (k=0; k<dim; k++) {
        B[j][k] *= mag2;
        }
      }
      for (k=0; k<dim; k++) {
        B[i][k] *= mag1;
      } 
    }
	return true;
}

//compute GramSchmidt orthogonalisation without normalisation
void GramSchmidt(int dim, int start, double **B, double *Mu) {
  int i, j, k; 
  //double mu_ij;
  double vec1[dim]; //store values to subtract from initial vectors
  double mag1;
  double mag2;
  bool orth_check = false;
  while (!orth_check) {
    orth_check = true;
    //iterate through the initial vectors
    for (i=fmax(start, 1); i<dim; i++) { 
      mag1 = sqrt(InnerProduct(dim, B[i], B[i]));
      for (j=0; j<dim; j++) {
        vec1[j] = 0;
        B[i][j] /= mag1; //normalising vectors before computing inner proucts helps to reduce inaccuracies caused by double calculations
      }
      //iterate through the previous vectors
      for (j=0; j<i; j++) { 
        mag2 = sqrt(InnerProduct(dim, B[j], B[j]));
        for (k=0;k<dim;k++) {
          B[j][k] /= mag2; //normalise before inner product
        }
        //if (j<2 && i<3) {
          //printf("j: %d\t", j);
          //printf("ip check3: %.40f\n", InnerProduct(dim, B[i], B[j]) / InnerProduct(dim, B[j], B[j]) - B[i][0]/B[j][0]);
        
          //for (k=0;k<5;k++) {
            //printf("%.4f\t", B[i][k]);
            //printf("%.4f\n", B[j][k]);
          //}
        //}
        Mu[(i-1)*i/2+j] = InnerProduct(dim, B[i], B[j]);///InnerProduct(dim, B[j], B[j]);
        for (k=0; k<dim; k++) {
          vec1[k] += Mu[(i-1)*i/2+j] * B[j][k] * mag1; //add the dot_product times the jth normalised vector 
          
          B[j][k]*=mag2;
        }
        
        //if (j==0 && i<5) {
        //  printf("mu_ij: %.40f\n", mu_ij);
        //}
        
        
      }
      //subtract from the ith initial vector
      for (k=0; k<dim; k++) {
        B[i][k]*=mag1; //reset B[i] to unnormalised version
        B[i][k] -= vec1[k];
      }
    }
	  //orth_check = CheckOrth(dim, start, B);
  }
}
  
//when A gets updated, recompute B to be the GramSchmidt orthogonalised version of the updated A
void update_matrices(int dim, int start, double **A, double **B, double *Mu) {
  int i, j;

  //set B to equal A
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      B[i][j] = A[i][j];  
    }
  }
  
  GramSchmidt(dim, start, B, Mu);
  //printf("Yes\n");
  for (i=0; i<dim; i++) {
    for (j=0; j<i; j++) {
      printf("%.4f\t", InnerProduct(dim, B[i], B[j]));
    }
  }
  printf("\n\n");
}


/// Lenstra–Lenstra–Lovász reduce the input matrix A
void LLL(double delta, int dim, double **A, double **B, double *Mu) {
  
  int i, j, k; //initialise variables i, j, k
  //double B[dim][dim];
  //set be to be equal to A
  //for (i=0; i<dim; i++) { 
  //  for (j=0; j<dim; j++) { 
  //    B[i][j] = A[i][j]; 
  //    } 
  //}

  //GramSchmidt orthogonslise B
  //GramSchmidt(dim, 0, B); 
  update_matrices(dim, 0, A, B, Mu);
  //printf("B: \n");
  //for (i=0; i<dim; i++) {
    //for (j=0; j<dim; j++) {
    //  printf("%.4f\t", B[i][j]);
    //}
    //printf("\n");
  //}
  
  k = 1;
  //double mu_kj;
  //double mu_k_kminus1; 
  bool zero_check = false; //checks if any vectors are reduced to 0, this can happen if the input vectors are linearly dependent and can result in an infinite loop
  //iterate through the LLL Reduction steps until:
  //(B[k] . B[k]) > (delta - mu_k_k-1) * (B[k-1] . B[k-1]) for every k, and
  //mu_kj<=0.5 for all k, j<k
  int m = 0;
  while (k<dim) {
    //reduce the kth vector until for all j<k, mu_kj<=0.5
    for (j=k-1; j>=0; j--) {
      //printf("IP: %.4f\n", InnerProduct(dim, 
      Mu[(k-1)*k/2+j] = InnerProduct(dim, A[k], B[j])/InnerProduct(dim, B[j], B[j]); 
      if (fabs(Mu[(k-1)*k/2+j]) > 0.5) {
        zero_check = true;
        for (i=0; i<dim; i++) {
          A[k][i] -= round(Mu[(k-1)*k/2+j]) * A[j][i];
          if (A[k][i] != 0) zero_check = false;
        }
        if (zero_check == true) {
          printf("Error: input vectors are linearly dependent\n");
          exit(1);
        }
        update_matrices(dim, 0, A, B, Mu);
      }
    }
    //LLL basis reduction requires (B[k] . B[k]) > (delta - mu_k_k-1) * (B[k-1] . B[k-1]) for every k
    Mu[(k-1)*k/2+k-1] = InnerProduct(dim, A[k], B[k-1])/InnerProduct(dim, B[k-1], B[k-1]); 
    if (InnerProduct(dim, B[k], B[k]) > ((delta - (Mu[(k-1)*k/2+k-1]*Mu[(k-1)*k/2+k-1])) * InnerProduct(dim, B[k-1], B[k-1]))) {
      k+=1;
        }
    else {
      //swap A[k] and A[k-1]
      for (i=0; i<dim; i++) {
        A[k][i] += A[k-1][i];
        A[k-1][i] = A[k][i] - A[k-1][i];
        A[k][i] -= A[k-1][i];
      }
      update_matrices(dim, 0, A, B, Mu); 
      k = fmax(k-1, 1);          
    }
    m++;
    //printf("%d\n", m);
    if (m % 100000 == 0) { //need to improve this
      printf("While loop failed\n");
        exit(1);
    }
  }
  printf("m: %d\n", m);
  //for (i=0;i<dim;i++) {
	//	for (j=0;j<dim;j++) {
	//		printf("%.4f\t", B[i][j]);
	//	}
	//	printf("\n");
	//}
	//printf("Address of B after LLL function: %p\n", (void *)B);
	//printf("dim: %d\n", dim);
}
  
  
  


