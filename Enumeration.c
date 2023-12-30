#include "LLL_Reduction.h"
#include <math.h>
#include <stdio.h>


//double InnerProduct(int dim, int *arr1, int *arr2) {
 // double sum1;
  //int i;
  //for (i=0; i<dim; i++) {
   // sum1 += arr1[i]*arr2[i];
  //}
  //return sum1;
//}

double ShortestVector(int dim, double (*A)[dim]) {
	int i, j, k;
	double shortest_vector = sqrt(InnerProduct(dim, A[0], A[0])); //keeps track of current shortest vector
	double current_norm; //stores current norm we're calculating below, may be good to release this memory after following for loop
	for (i=1; i<dim; i++) {
		current_norm = sqrt(InnerProduct(dim, A[i], A[i])); //calculates norm of each vector in A, then compares o shortest vector
		if (current_norm < shortest_vector) {
			shortest_vector = current_norm;
		}
	}
	printf("shortest basis vector: %.4f\n", shortest_vector);
	double Mu[dim][dim]; //stores the mu_i_j values
	double GS_norms[dim]; //stores the norm of each GS vector
	double sum1[dim]; //may be able to release this memory after following for loop
	for (i=1; i<dim; i++) { //change A to GS vectors, calculate Mu_i_j values, and calculate GS_norms
		GS_norms[i-1] = InnerProduct(dim, A[i-1], A[i-1]);
		//printf("GS_norms: %.4f\n", GS_norms[i-1]);
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
	GS_norms[dim-1] = InnerProduct(dim, A[dim], A[dim]);
	int x[dim];
	double l[dim];
	for (j=0; j<dim; j++) {
		x[j] = l[j] = 0;
	}
	double sum2;
	double sum3;
	i=0;
	int m = 0;
	while (i<dim) {
		//printf("i: %d\n", i);
		//for (j=0; j<dim; j++) {
			//printf("x: %d\t", x[j]);
		//}
		//printf("\n");
		sum3 = sum2 = 0;
		for (j=dim-1; j>=i; j--) {
			if (j<dim-1) {
				sum2 += x[j+1] * Mu[j+1][j];
			}
			//sum2 = 0;
			//for (k=j+1; k<dim; k++) {
			//	sum2 += x[k] * Mu[k][j];
			//}
			//printf("x[j] + sum2: %.4f\n", x[j] + sum2);
			l[j] = (x[j] + sum2) * (x[j] + sum2) * GS_norms[j];
			if (j==i) {
				printf("GS_norms[i]: %.4f\n", GS_norms[i]);
				printf("sum2: %.4f\n", sum2);
				printf("sum3: %.4f\n", sum3);
				printf("xi...: %.4f\n", (x[i] + sum2)*sqrt(GS_norms[i]) - sqrt(shortest_vector - sum3));
			}
			
		}
		for (j=0; j<dim; j++) {
			printf("l: %.4f\t", l[j]);
		}
		printf("\n");
		for (j=0; j<dim; j++) {
			sum3 += l[j];
		}
		if (i==0 && sum3 < shortest_vector) {
			if (sum3 != 0) {
				shortest_vector = sum3;
			}
			x[0] += 1;
		}
		else if (i!=0) {
			sum3 = 0;
			for (j=i; j<dim; j++) {
				sum3 += l[j];
			}
			//printf("sum3: %.4f\n", sum3);
			if (sum3<shortest_vector) {
				printf("i, sum3, x[i]: %d %.4f %d\n", i, sum3, x[i]);
				i -= 1;
				x[i] = ceil(- sum2 - x[i+1]*Mu[i+1][i] - sqrt((shortest_vector - sum3)/GS_norms[i]));
				printf("i2, sum2, x[i]: %d %.4f %d\n", i, sum2 + x[i+1]*Mu[i+1][i], x[i]);
			}
		}
		if (sum3 >= shortest_vector) {
			i += 1;
			x[i] += 1;
		}
		m+=1;
		if (m==100) {
			break;
		}
	}
	return shortest_vector;
}
