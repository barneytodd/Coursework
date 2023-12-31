#include "LLL_Reduction.h"
#include <math.h>
#include <stdio.h>


//double InnerProduct(int dim, double *arr1, double *arr2) {
//  double sum1;
//  int i;
 // for (i=0; i<dim; i++) {
  //  sum1 += arr1[i]*arr2[i];
  //}
  //return sum1;
//}

double ShortestVector(int dim, double (*A)[dim]) {
	int i, j, k;
	double shortest_vector = sqrt(InnerProduct(dim, A[0], A[0])); //keeps track of current shortest vector
	double current_norm; //stores current norm we're calculating below, may be good to release this memory after following for loop
	for (i=1; i<dim; i++) {
		current_norm = sqrt(InnerProduct(dim, A[i], A[i])); //calculates norm of each vector in A, then compares to shortest vector
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

	//for (i=0;i<dim;i++) {
	//	for (j=i+1;j<dim;j++) {
	//		printf("innerproduct, i, j: %.4f %d %d\n", InnerProduct(dim, A[i], A[j]), i, j);
	//	}
	//}
	
	GS_norms[dim-1] = InnerProduct(dim, A[dim], A[dim]);
	printf("GS: %.4f\n", GS_norms[9]);
	int x[dim]; //x counts how many of each basis vector we're using
	double l[dim]; //l counts the total contribution of all the used vectors in the direction of each GS vector, squared
	for (j=0; j<dim; j++) {
		x[j] = l[j] = 0;
	}
	double sum2; //sum2 is the sum(j>i) of x[j] * Mu[j][i]
	double sum3; //sum3 is the sum of l[j]'s, also equal to the squared norm of the current vector
	i=0;
	//int m = 0;
	
	while (i<dim) { //begin the enumeration loop
		//printf("i: %d\n", i);
		//for (j=0; j<dim; j++) {
			//printf("x: %d\t", x[j]);
		//}
		//printf("\n");
		sum2 = 0;
		for (j=dim-1; j>=i; j--) { //calculate the l[j] values from i upwards
			//if (j<dim-1) { //dynamically calculate sum2
			//	sum2 += x[j+1] * Mu[j+1][j];
			//}
			sum2 = 0;
			for (k=j+1; k<dim; k++) {
				sum2 += x[k] * Mu[k][j];
			}
			//printf("x[j] + sum2: %.4f\n", x[j] + sum2);
			l[j] = (x[j] + sum2) * (x[j] + sum2) * GS_norms[j]; //might need to include all j not just j>i
			if (j==i) {
				//printf("l[j]...: %.4f\n", l[j] - (shortest_vector-sum3));
				//printf("GS_norms[i]: %.4f\n", GS_norms[i]);
				//printf("sum2: %.4f\n", sum2);
				//printf("sum3: %.4f\n", sum3);
				//printf("x1...: %.4f\n", ((- sum2 - sqrt((shortest_vector - sum3)/GS_norms[i]))+sum2)*sqrt(GS_norms[i])- sqrt(shortest_vector - sum3));
				//printf("x2...: %.4f\n", ((- sum2 - sqrt((shortest_vector - sum3)/GS_norms[i]))+sum2)-sqrt((shortest_vector - sum3)/GS_norms[i]));
				//printf("x3...: %.4f\n", ((- sum2 - sqrt((shortest_vector - sum3)/GS_norms[i]))+sum2)*sqrt(GS_norms[i])-sqrt((shortest_vector - sum3)));
				//printf("xi...: %.4f\n", (x[i] + sum2)*sqrt(GS_norms[i]) - sqrt(shortest_vector - sum3));
			}
			
		}
		//for (j=0; j<dim; j++) {
		//	printf("l: %.4f\t", l[j]);
		//}
		//printf("\n");
		sum3 = 0;
		for (j=i; j<dim; j++) {
			sum3 += l[j];
		}
		if (sum3 < shortest_vector*shortest_vector) {
			if (i==0) {
				if (sum3 != 0) {
					shortest_vector = sqrt(sum3);
					printf("shortest_vector: %.4f\n", shortest_vector);
					for (j=0; j<dim; j++) {
						printf("%d\t", x[j]);
					}
					printf("\n");
					//printf("x[9]: %
					printf("max x[9]: %.4f\n", shortest_vector*shortest_vector/GS_norms[dim-1]);
				}
				x[0] += 1;
			}
			else {
				//we want to set x[i] to the minimum integer such that l[i] < A-sum3


				
				//printf("i, sum3, x[i]: %d %.4f %d\n", i, sum3, x[i]);
				i -= 1;
				//x[i] = ceil(- sum2 - x[i+1]*Mu[i+1][i] - sqrt((shortest_vector - sum3)/GS_norms[i]));
				x[i] = round(- sum2 - x[i+1]*Mu[i+1][i]);
				//printf("i2, sum2, x[i]: %d %.4f %d\n", i, sum2 + x[i+1]*Mu[i+1][i], x[i]);
			}
		}
		else {

			//we might want to reset all the x[i]'s to 0 when x[9] gets updated, also not sure about max x[9]

			
			i += 1;
			x[i] += 1;
		}
		//if (i==dim-1) {
		//	printf("x[9]: %d\n", x[dim-1]);
		//	
		//}
		//m+=1;
		//if (m==100) {
		//	break;
		//}
	}
	return shortest_vector;
}
