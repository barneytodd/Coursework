#include "LLL_Reduction.h"


double InnerProduct(int dim, int *arr1, int *arr2) {
  double sum1;
  int i;
  for (i=0; i<dim; i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sqrt(sum1);
}

double ShortestVector(int dim, double (*A)[dim]) {
	int i, j, k;
	double shortest_vector = InnerProduct(dim, A[0], A[0]); //keeps track of current shortest vector
	double current_norm; //stores current norm we're calculating below, may be good to release this memory after following for loop
	for (i=1; i<dim; i++) {
		current_norm = InnerProduct(dim, A[i], A[i]); //calculates norm of each vector in A, then compares o shortest vector
		if (current_norm < shortest_vector) {
			shortest_vector = current_norm;
		}
	}
	double Mu[dim][dim]; //stores the mu_i_j values
	double GS_norms[dim]; //stores the norm of each GS vector
	double sum1[dim]; //may be able to release this memory after following for loop
	for (i=1; i<dim; i++) { //change A to GS vectors, calculate Mu_i_j values, and calculate GS_norms
		GS_norms[i-1] = InnerProduct(dim, A[i-1], A[i-1]);
		for (j=0; j<dim; j++) {
			sum1[j] = 0;
		}
		for (j=0; j<i; j++) {
			Mu[i][j] = InnerProduct(dim, A[i], A[j])/GS_norms[j];
			for (k=0; k<dim; k++) {
				sum1[k] += Mu[i][j] * A[j][k];
			}
		}
		for (k=0; k<dim; k++) {
			A[i][k] -= sum1[k];
		}	
	}
	int x[dim];
	double l[dim];

	double sum2;
	i=0;
	while (i<dim) {
		sum2 = 0;
		for (j=i+1; j<dim; j++) {
			sum2 += x[j] * Mu[j][i];
		}
		l[i] = (x[i] + sum2) * (x[i] + sum2) * GS[i];
		sum2 = 0;
		for (j=0; j<dim; j++) {
			sum2 += l[j];
		}
		if (i==0 && sum2 < shortest_vector) {
			shortest_vector = sum2;
			x[0] += 1;
		}
		else {
			sum2 = 0;
			for (j=i; j<dim; j++) {
				sum2 += l[j];
			}
			if (sum2<shortest_vector) {
				i -= 1;
				
			}
		}
	}
}
