#include "LLL_Reduction.h"
#include <math.h>
#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

struct ThreadArgs {
  int num;
  int dim;
  double *GS_norms;
  double *Mu;
  double *shortest_vector;
  int *max_num;
	pthread_mutex_t *lock;
};

void *Enumerate(void *args) {
  //need num for x[dim-1], GS_norms, Mu, dim, pointer to shortest_vector, pointer to max x[dim-1], create own x, create own l
  struct ThreadArgs *thread_args = (struct ThreadArgs *)args;

  int x[thread_args->dim], i, j, k;
  double l[thread_args->dim];

  for (i=0; i<thread_args->dim-1; i++) {
    x[i] = 0;
  }
  x[thread_args->dim-1] = thread_args->num;
  i = thread_args->dim-1;

  double sum2; //stores the sum of x[j] * Mu[j][i] for j>i
  double sum3;
  int m = 0;
  while (*thread_args->max_num > thread_args->num) { 
		sum2 = 0;
		//calculate the l[j] values from i upwards
		for (j=thread_args->dim-1; j>=i; j--) { 
			sum2 = 0;
			//sum the contribution of each vector in the direction of the ith GS vector 
			for (k=j+1; k<thread_args->dim; k++) {
				sum2 += x[k] * thread_args->Mu[(k-1)*k/2+j]; //Mu[k][j]
			}
			l[j] = (x[j] + sum2) * (x[j] + sum2) * thread_args->GS_norms[j]; 	
		}

		//sum the l[j] values for j>=i
		sum3 = 0;
		for (j=i; j<thread_args->dim; j++) {
			sum3 += l[j];
		}
		
		if (sum3 < (*(thread_args->shortest_vector))*(*(thread_args->shortest_vector))) {
			//if i=0 and sum3 < (current shortest vector length)^2, we have a new shortest vector
			if (i==0) {
				if (sum3 != 0) {
					pthread_mutex_lock(thread_args->lock);
					*(thread_args->shortest_vector) = sqrt(sum3);
					printf("shortest_vector: %.4f\n", *(thread_args->shortest_vector));
					for (j=0; j<thread_args->dim; j++) {
						printf("%d\t", x[j]);
					}
					printf("\n");
					*(thread_args->max_num) = floor(*(thread_args->shortest_vector)/pow(thread_args->GS_norms[thread_args->dim-1], 0.5));
					printf("max x[dim-1]: %d\n", *(thread_args->max_num));
					pthread_mutex_unlock(thread_args->lock);
				}
				x[0] += 1;
			}
			//if i != 1, subtract 1 from i and then set the x[i] to be the minimum integer such that l[i] < shortest_vector^2 - sum3
			//if there is no such integer, add the 1 back to i and then add 1 to x[i]
			else {
				i -= 1;
				
				sum2 = 0;
				for (k=i+1; k<thread_args->dim; k++) {
					sum2 += x[k] * thread_args->Mu[(k-1)*k/2+i]; //Mu[k][i]
				}
				x[i] = round(- sum2); //the integer which minimises l[i]
				l[i] = ((double)x[i] + sum2) * ((double)x[i] + sum2) * thread_args->GS_norms[i]; 
				
				if (l[i] < (*(thread_args->shortest_vector)) * (*(thread_args->shortest_vector)) - sum3) {
					//subtract 1 from x[i] until l[i] is no longer < shortest_vector^2 - sum3
					//then add 1 to x[i] to make x[i] the minimum possible integer such that l[i] < shortest_vector^2
					do {
						
						x[i] -= 1;
						sum2 = 0;
						for (k=i+1; k<thread_args->dim; k++) {
							sum2 += x[k] * thread_args->Mu[(k-1)*k/2+i]; //Mu[k][i]
						}
						l[i] = (x[i] + sum2) * (x[i] + sum2) * thread_args->GS_norms[i]; 
						
						
					} while (l[i] < (*(thread_args->shortest_vector)) * (*(thread_args->shortest_vector)) - sum3);
					x[i] += 1; 
				}
				else {
					i+=1;
					x[i] += 1;
				}
			}
		}
		//if sum3 > shortest_vector^2, increase i by 1 and then increase x[i] by 1
		else {
			i += 1;
			x[i] += 1;
			if (i==thread_args->dim-1) {
				break;
			}
		}
	  m++;
	  if (m==1000000000000000) { //need to improve this
		  printf("thread: %d, infinite while loop\n", thread_args->num);
			exit(1);
	  }
  }
  pthread_exit(NULL);
}

//Enumerate the lattice to find the shortest vector
double ShortestVector1(int dim, double **A) {
	for (int i = 0; i < dim; ++i) {
		if (A[i] == NULL) {
		    perror("Input matrix does not have the correct dimensions");
		    exit(1);
		}
	}
	
	
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
	for (i=1; i<dim; i++) { 
		GS_norms[i-1] = InnerProduct(dim, A[i-1], A[i-1]);
		
		for (j=0; j<dim; j++) {
			sum1[j] = 0;
		}
		for (j=0; j<i; j++) {
			Mu[(i-1)*i/2+j] = InnerProduct(dim, A[i], A[j])/GS_norms[j];
			for (k=0; k<dim; k++) {
				sum1[k] += Mu[(i-1)*i/2+j] * A[j][k]; //equivalent to Mu[i][j] if Mu was a dim x dim array
			}
		}
		for (k=0; k<dim; k++) {
			A[i][k] -= sum1[k];
		}	
	}
	GS_norms[dim-1] = InnerProduct(dim, A[dim-1], A[dim-1]);
	for (i=0;i<dim;i++) {
		printf("%.4f\t", GS_norms[i]);
	}
	printf("\n");
	printf("shortest basis vector: %.4f\n", shortest_vector); 
	//int x[dim]; //counts how many of each basis vector we're using
	//double l[dim]; //stores the total contribution of all the used vectors in the direction of each GS vector, squared

	//initialise the x and l values to 0
	//for (j=0; j<dim; j++) {
	//	x[j] = l[j] = 0;
	//}
	
	//double sum2; //stores the sum of x[j] * Mu[j][i] for j>i
	//double sum3; //sum of l[j]'s, (the sum of all the l[j]'s from 1 to dim is equal to the squared norm of the current combination of basis vectors)
	//i=dim-1; //start with the last vector
	int max_num = floor(shortest_vector/pow(GS_norms[dim-1], 0.5));
	pthread_t threads[max_num];

	pthread_mutex_t lock;
	if (pthread_mutex_init(&lock, NULL) != 0) {
        	printf("Error: Mutex initialization failed\n");
        	exit(1);
    	}
	
	struct ThreadArgs args[max_num+1];
	int count = 0;
	for (i=0; i<=max_num; i++) {
		args[i].num = i;
		args[i].dim = dim;
		args[i].GS_norms = GS_norms;
		args[i].Mu = Mu;
		args[i].shortest_vector = &shortest_vector;
		args[i].max_num = &max_num;
		args[i].lock = &lock;
		if (pthread_create(&threads[i], NULL, &Enumerate, (void *)&args[i]) != 0) {
			printf("Error creating thread\n");
			exit(1);
		}
		printf("Thread: %d\n", i);
		count++;
	}
	for (i=0; i<count; i++) {
		pthread_join(threads[i], NULL);
	}
	pthread_mutex_destroy(&lock);
	return shortest_vector;
}
