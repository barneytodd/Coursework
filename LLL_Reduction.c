#include "LLL_Reduction.h"
#include "GeneralFunctions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

//compute GramSchmidt orthogonalisation without normalisation
void GramSchmidt(int dim, int start, double **B, double *Mu) {
  int i, j, k; 
  //double vec1[dim]; //store values to subtract from initial vectors
  double mag1; //stores the norm of each i vector in the following loop
  double mag2; //stores the norm of each j vector
	
	//iterate through the initial vectors
	for (i=fmax(start, 1); i<dim; i++) { 
		mag1 = sqrt(InnerProduct(dim, B[i], B[i]));
		for (k=0; k<dim; k++) {
			//vec1[k] = 0;
			B[i][k] /= mag1; //normalising vectors before computing inner proucts helps to reduce inaccuracies caused by double calculations
		}
		for (j=0; j<i; j++) {
			mag2 = sqrt(InnerProduct(dim, B[j], B[j]));
			for (k=0;k<dim;k++) {
		  	B[j][k] /= mag2; //normalise before inner product
			}
			Mu[(i-1)*i/2+j] = InnerProduct(dim, B[i], B[j])*mag1;
		}
		for (k=0; k<dim; k++) {
			B[i][k] *= mag1;
		}
		//iterate through the previous vectors
		for (j=0; j<i; j++) { 
			mag2 = sqrt(InnerProduct(dim, B[j], B[j]));
			//for (k=0;k<dim;k++) {
		  	//B[j][k] /= mag2; //normalise before inner product
			//}
			//Mu[(i-1)*i/2+j] = InnerProduct(dim, B[i], B[j])*mag1;///InnerProduct(dim, B[j], B[j]);
			for (k=0; k<dim; k++) {
		  		//vec1[k] += Mu[(i-1)*i/2+j] * B[j][k]; //add the dot_product times the jth normalised vector 
				if (k==0) {
				printf("Mu: %.4f, B[i]: %.4f, B[j]: %.4f\n", Mu[(i-1)*i/2+j], B[i][k], B[j][k]);
				}
				B[i][k] -= Mu[(i-1)*i/2+j] * B[j][k];
				B[j][k] *= mag2; //reset B[j] to original values
			} 
			Mu[(i-1)*i/2+j] /= mag2; //reset Mu[i][j] to normal value     
		}
				
		//subtract from the ith initial vector
		//for (k=0; k<dim; k++) {
			//B[i][k]*=mag1; //reset B[i] to original version
			//B[i][k] -= vec1[k];
		//}
	}
	printf("B\n");
	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			printf("%.4f ", B[i][j]);
		}
		printf("\n");
	}
}
  
//when A gets updated, recompute B to be the GramSchmidt orthogonalised version of the updated A, along with the associated Mu values
void update_matrices(int dim, int start, double **A, double **B, double *Mu) {
  int i, j;
  //set B to equal A for the vectors after the one that has just changed
  for (i=start; i<dim; i++) {
    for (j=0; j<dim; j++) {
      B[i][j] = A[i][j];  
    }
  }
  GramSchmidt(dim, start, B, Mu);
}


/// Lenstra–Lenstra–Lovász reduce the input matrix A
void LLL(double delta, int dim, double **A, double **B, double *Mu) {

  int i, j, k; //initialise variables i, j, k

	//check the size of A and B
	for (i = 0; i < dim; i++) {
		if (A[i] == NULL || B[i] == NULL) {
			printf("Error: Input matrices for ShortestVector do not have the correct dimensions");
			FreeMatrix(i, &A);
			FreeMatrix(i, &B);
			free(Mu);
			Mu = NULL;		
			exit(1);
		}
	}
	
  update_matrices(dim, 0, A, B, Mu); 
    
  k = 1;
  bool zero_check = false; //checks if any vectors are reduced to 0, this can happen if the input vectors are linearly dependent and can result in an infinite loop
  int m = 0;
	
	//iterate through the LLL Reduction steps until:
  //(B[k] . B[k]) > (delta - mu_k_k-1) * (B[k-1] . B[k-1]) for every k, and
  //mu_kj<=0.5 for all k, j<k
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
					FreeMatrix(dim, &A);
					FreeMatrix(dim, &B);
					free(Mu);
					Mu = NULL;		
          exit(1);
        }
        update_matrices(dim, k, A, B, Mu);
				m=0;
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
        B[k][i] = A[k-1][i]; //B[k] is about to get updated, so we can use this as a temporary variable
        A[k-1][i] = A[k][i];
        A[k][i] = B[k][i];			
      }
      update_matrices(dim, k-1, A, B, Mu); 
      k = fmax(k-1, 1);  
			m++;
    }

		//(dim-1)*dim/2 is the number of swaps required to completely reverse the list of vectors
		//therefore any more than this without reducing any of the vectors means an error has occured
    if (m > (dim-1)*dim/2) { 
      printf("Error: LLL while loop failed\n");
			FreeMatrix(dim, &A);
			FreeMatrix(dim, &B);
			free(Mu);
			Mu = NULL;		
    	exit(1);
    }
  }
}
  
  
  


