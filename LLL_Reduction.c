#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "LLL_Reduction.h"
#include <stdbool.h>
#include <limits.h>

double Multiply(double *num1, double *num2) {
	printf("%.4f\n", *num1);
	if (*num1<pow(10, 15) && *num2<pow(10, 15)) {
		return (*num1)*(*num2);
	}
	else {
		int i, j;
		char *str1 = (char *)malloc(16*sizeof(char));
		char *str2 = (char *)malloc(1*sizeof(char));
		int len1 = snprintf(str1, 16, "%.0f", fmax(*num1, *num2));
		int len2 = snprintf(str1, 4, "%.0f", fmin(*num1, *num2));
		if (len1 > 16) {
			str1 = (char *)realloc(str1, (len1+1) * sizeof(char));
        		snprintf(str1, len1+1, "%.0f", fmax(*num1, *num2)); 
		}
		if (len2 > 2) {
			str2 = (char *)realloc(str2, (len2+1) * sizeof(char));
        		snprintf(str2, len2+1, "%.0f", fmin(*num1, *num2)); 
		}
		char *result = (char *)malloc((len1+len2+2)*sizeof(char));
		for (i=0; i<len1+len2+1; i++) {
			result[i] = '0';
		}
		printf("str1: %s\t str2: %s\n", str1, str2);
		int carry;
		int current;
		for (i=0; i<len2; i++) {
			for (j=0; j<len1; j++) {
				carry = 0;
				current = (int)str2[i]*(int)str1[j]+result[i+j+1]+carry;
				while (current >= 10) {
					current -= 10;
					carry ++;
				}
				result[i+j+1] = current + '0';
				//snprintf(&result[i+j+1], 1, "%d", current);
			}
		}
		return strtod(result, NULL);
	}
	
}



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
        if (fabs(InnerProduct(dim, B[i], B[j])) > 0.0001) {
          //printf("Failed: %.4f\n", fabs(InnerProduct(dim, B[i], B[j])));
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
	  orth_check = CheckOrth(dim, start, B);
  }
}
  
//when A gets updated, recompute B to be the GramSchmidt orthogonalised version of the updated A
void update_matrices(int dim, int start, double **A, double **B, double *Mu) {
  int i, j;

  //set B to equal A for the vectors after the one that has just changed
  for (i=start; i<dim; i++) {
    for (j=0; j<dim; j++) {
      B[i][j] = A[i][j];  
    }
  }
  
  GramSchmidt(dim, start, B, Mu);
  //printf("Yes\n");
  //for (i=0; i<dim; i++) {
    //for (j=0; j<i; j++) {
      //printf("%.4f\t", InnerProduct(dim, B[i], B[j]));
    //}
  //}
  //printf("\n\n");
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
	double array[dim][dim];
	for (i=0;i<dim;i++) {
		for (j=0; j<dim; j++) {
			if (i==j) {
				array[i][i] = 1;
			}
			else {
				array[i][j] = 0;
			}
		}
	}
  while (k<dim) {
    //reduce the kth vector until for all j<k, mu_kj<=0.5
    for (j=k-1; j>=0; j--) {
      //printf("IP: %.4f\n", InnerProduct(dim, 
      Mu[(k-1)*k/2+j] = InnerProduct(dim, A[k], B[j])/InnerProduct(dim, B[j], B[j]); 
      if (fabs(Mu[(k-1)*k/2+j]) > 0.5) {
        zero_check = true;
        for (i=0; i<dim; i++) {
		
		printf("%.4f\n", Multiply(round(Mu[(k-1)*k/2+j]), A[j][i]);
		exit(1);
					printf("k: %d, j: %d, A[k][i]: %.4f, round(Mu): %.1f, A[j][i]: %.4f, Mu*AJ: %.4f\n", k, j, A[k][i], round(Mu[(k-1)*k/2+j]), A[j][i], round(Mu[(k-1)*k/2+j]) * A[j][i]);
					printf("array[k][i]: %.4f, array[j][i]: %.4f\n", array[k][i], array[j][i]);
          A[k][i] -= round(Mu[(k-1)*k/2+j]) * A[j][i];
					printf("A[k][i]: %.4f\n", A[k][i]);
					array[k][i] -= round(Mu[(k-1)*k/2+j]) * array[j][i];
					printf("array[k][i]: %.4f\n", array[k][i]);
          if (A[k][i] != 0) zero_check = false;
        }
        if (zero_check == true) {
          printf("Error: input vectors are linearly dependent\n");
          exit(1);
        }
        update_matrices(dim, k, A, B, Mu);
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

				B[k][i] = array[k-1][i];
				array[k-1][i] = array[k][i];
				array[k][i] = B[k][i];
      }
      update_matrices(dim, k-1, A, B, Mu); 
      k = fmax(k-1, 1);          
    }
		
    m++;
    //printf("%d\n", m);
		if (m % 10 ==0) {
			printf("A\n");
	    for (i=0;i<1;i++) {
				for (j=0;j<dim;j++) {
					printf("%.4f ", A[i][j]);
				}
				printf("\n");
			}
			printf("array\n");
	    for (i=0;i<1;i++) {
				for (j=0;j<dim;j++) {
					printf("%.4f ", array[i][j]);
				}
				printf("\n");
			}
		}
    if (m % 30 == 0) { //need to improve this
      printf("While loop failed\n");
    	exit(1);
    }
  }
	printf("shortest\n");
	for (i=0; i<dim; i++) {
		printf("%f, ", array[0][i]);
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
  
  
  


