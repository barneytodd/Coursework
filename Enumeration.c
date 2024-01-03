#include "LLL_Reduction.h"
#include <math.h>
#include <stdio.h>


//Enumerate the lattice to find the shortest vector
double ShortestVector(int dim, double **A) {
	int i, j, k;
	double shortest_vector = sqrt(InnerProduct(dim, A[0], A[0])); //keeps track of current shortest vector
	double current_norm; //stores the norm of each basis vector

	//find the shortest basis vector and set shortes_vector to be equal to its norm
	for (i=1; i<dim; i++) {
		current_norm = sqrt(InnerProduct(dim, A[i], A[i])); //calculates norm of each vector in A, then compares to shortest vector
		if (current_norm < shortest_vector) {
			shortest_vector = current_norm;
		}
	}
	
	double Mu[(dim-1)*dim/2]; //stores the mu_i_j values in a lower triangular matrix
	double GS_norms[dim]; //stores the norm of each GramSchidt orthogonalised vector
	
	//GramSchmidt orthogonalise A, and store the mu values and GS norms
	double sum1[dim]; 
	for (i=1; i<dim; i++) { //change A to GS vectors, calculate Mu_i_j values, and calculate GS_norms
		GS_norms[i-1] = InnerProduct(dim, A[i-1], A[i-1]);
		
		for (j=0; j<dim; j++) {
			sum1[j] = 0;
		}
		for (j=0; j<i; j++) {
			Mu[(i-1)*i/2+j] = InnerProduct(dim, A[i], A[j])/GS_norms[j];
			for (k=0; k<dim; k++) {
				sum1[k] += Mu[(i-1)*i/2+j] * A[j][k]; //equivalent to Mu[i][j]
			}
		}
		for (k=0; k<dim; k++) {
			A[i][k] -= sum1[k];
		}	
	}
	GS_norms[dim-1] = InnerProduct(dim, A[dim-1], A[dim-1]);
	
	int x[dim]; //counts how many of each basis vector we're using
	double l[dim]; //stores the total contribution of all the used vectors in the direction of each GS vector, squared

	//initialise the x and l values to 0
	for (j=0; j<dim; j++) {
		x[j] = l[j] = 0;
	}
	
	double sum2; //stores the sum of x[j] * Mu[j][i] for j>i
	double sum3; //sum of l[j]'s, also equal to the squared norm of the current combination of basis vectors
	i=dim-1;
	
	//enumeration the lattice
	while (i<dim) { 
		sum2 = 0;
		//calculate the l[j] values from i upwards
		for (j=dim-1; j>=i; j--) { 
			sum2 = 0;
			//sum the contribution of each vector in the direction of the ith GS vector 
			for (k=j+1; k<dim; k++) {
				sum2 += x[k] * Mu[(k-1)*k/2+j]; //Mu[k][j]
			}
			l[j] = (x[j] + sum2) * (x[j] + sum2) * GS_norms[j]; 	
		}

		//sum the l[j] values for j>=i
		sum3 = 0;
		for (j=i; j<dim; j++) {
			sum3 += l[j];
		}
		
		if (sum3 < shortest_vector*shortest_vector) {
			//if i=0 and sum3< (current shortest vector length)^2, we have a new shortest vector
			if (i==0) {
				if (sum3 != 0) {
					shortest_vector = sqrt(sum3);
					printf("shortest_vector: %.4f\n", shortest_vector);
					for (j=0; j<dim; j++) {
						printf("%d\t", x[j]);
					}
					printf("\n");
					
					printf("max x[9]: %.4f\n", pow(shortest_vector*shortest_vector/GS_norms[dim-1], 0.5));
					
					printf("\n");
				}
				x[0] += 1;
			}
			else {
				//we want to set x[i] to the minimum integer such that l[i] < A-sum3


				
				i -= 1;
				
				sum2 = 0;
				for (k=i+1; k<dim; k++) {
					sum2 += x[k] * Mu[(k-1)*k/2+i]; //Mu[k][i]
				}
				x[i] = round(- sum2);
				
				//l[i] = ((double)x[i] + sum2) * ((double)x[i] + sum2) * GS_norms[i]; 
				
				do {
					
					x[i] -= 1;
					sum2 = 0;
					for (k=i+1; k<dim; k++) {
						sum2 += x[k] * Mu[(k-1)*k/2+i]; //Mu[k][i]
					}
					l[i] = (x[i] + sum2) * (x[i] + sum2) * GS_norms[i]; 
					
					
				} while (l[i] < shortest_vector * shortest_vector - sum3);
				x[i] += 1; 
			}
		}
		else {

			

			
			i += 1;
			x[i] += 1;
			if (i==dim-1) {
				printf("x[dim-1]: %d\n", x[dim-1]);
			}
		
		}
		
		//printf("\n");
	}
	return shortest_vector;
}
