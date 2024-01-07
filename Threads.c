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
	//restructure the arguments into the form of the struct above
  struct ThreadArgs *thread_args = (struct ThreadArgs *)args;

  int x[thread_args->dim], i, j, k; //x stores the number of each basis vector used to reach each lattice point
  double l[thread_args->dim]; //stores the total contribution squared, of the combination of basis vectors stored in x, in the direction of the ith GS vector

  for (i=0; i<thread_args->dim-1; i++) {
    x[i] = 0;
  }
  x[thread_args->dim-1] = thread_args->num;
  i = thread_args->dim-1;
	printf("num: %d\n", thread_args->num);
  double sum2; //stores the sum of x[j] * Mu[j][i] for j>i
  double sum3; //stores the sum (j>i) of the l[j] values
  int m = 0;
	int max_its = pow(2, thread_args->dim); //an estimate of the maximum number of iterations required
	
	//max_num may get updated by the other threads
	//in the case that max_num falls below num, we can exit this thread
  while (*thread_args->max_num > thread_args->num) { 
		//calculate the l[j] values from i upwards
		for (j=thread_args->dim-1; j>=i; j--) { 
			sum2 = 0;
			//sum the contribution of each vector in the direction of the ith GS vector 
			for (k=j+1; k<thread_args->dim; k++) {
				sum2 += x[k] * thread_args->Mu[(k-1)*k/2+j]; 
			}
			l[j] = (x[j] + sum2) * (x[j] + sum2) * thread_args->GS_norms[j]; 	
		}
		if (thread_args->num == 3) {
			printf("i: %d\n", i);
		}
		//sum the l[j] values for j>=i
		sum3 = 0;
		for (j=i; j<thread_args->dim; j++) {
			sum3 += l[j];
		}
	  	if (thread_args->num == 3) {
			printf("sum3: %.4f, %d\n", sum3, i);
			printf("l[8]: %.4f, l[9]: %.4f\n", l[8], l[9]);
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
				x[0]++;
			}
			//if i != 1, subtract 1 from i and then set x[i] to be the minimum integer such that l[i] < shortest_vector^2 - sum3
			//if there is no such integer, add the 1 back to i and then add 1 to x[i]
			else {
				i--;
				sum2 = 0;
				for (k=i+1; k<thread_args->dim; k++) {
					sum2 += x[k] * thread_args->Mu[(k-1)*k/2+i]; 
				}
				x[i] = round(- sum2); //the integer which minimises l[i]
				l[i] = ((double)x[i] + sum2) * ((double)x[i] + sum2) * thread_args->GS_norms[i]; 
				
				if (l[i] < (*(thread_args->shortest_vector)) * (*(thread_args->shortest_vector)) - sum3) {
					//subtract 1 from x[i] until l[i] is no longer < shortest_vector^2 - sum3
					//then add 1 to x[i] to make x[i] the minimum possible integer such that l[i] < shortest_vector^2 - sum3
					do {
						x[i]--;
						sum2 = 0;
						for (k=i+1; k<thread_args->dim; k++) {
							sum2 += x[k] * thread_args->Mu[(k-1)*k/2+i];
						}
						l[i] = (x[i] + sum2) * (x[i] + sum2) * thread_args->GS_norms[i]; 
					} while (l[i] < (*(thread_args->shortest_vector)) * (*(thread_args->shortest_vector)) - sum3);
					
					x[i]++; 
				}
				else {
					i+=1;
					x[i]++;
				}
			}
		}
		//if sum3 > shortest_vector^2, increase i by 1 and then increase x[i] by 1
		else {
			i++;
			x[i]++;
			//if we reach the point of incrementing x[dim-1], then we have tested all possibilities with x[dim-1] = num, therefore we can exit the thread
			if (i==thread_args->dim-1) {
				break;
			}
		}
		m++;
	  	if (thread_args->num == 3) {
		  	for (i=0; i<thread_args->dim; i++) {
				printf("%d ", x[i]);
			}
			printf("\n");
		
		  if (m > 100) {
				//need to think about this
		//	  printf("thread: %d, infinite while loop\n", thread_args->num);
				exit(1);
		  }
	  }
  }
  pthread_exit(NULL);
}

//Enumerate the lattice to find the shortest vector
double ShortestVector1(int dim, double **A, double **B, double *Mu) {
	int i;
	
	for (i = 0; i < dim; ++i) {
		if (A[i] == NULL) {
		    perror("Input matrix does not have the correct dimensions");
		    exit(1);
		}
	}
		
	double shortest_vector = sqrt(InnerProduct(dim, A[0], A[0])); //keeps track of current shortest vector
	double current_norm; //stores the norm of each basis vector

	//find the shortest basis vector and set shortes_vector to be equal to its norm
	for (i=1; i<dim; i++) {
		current_norm = sqrt(InnerProduct(dim, A[i], A[i])); //calculates norm of each vector in A, then compares to shortest vector
		if (current_norm < shortest_vector) {
			shortest_vector = current_norm;
		}
	}
	
	printf("shortest basis vector: %.4f\n", shortest_vector);
	
	double GS_norms[dim]; //stores the norm of each GramSchidt orthogonalised vector
	update_matrices(dim, 0, A, B, Mu);
	for (i=0;i<dim;i++) {
		GS_norms[i] = InnerProduct(dim, B[i], B[i]);
	}

	printf("GS\n");
	int max_num = floor(shortest_vector/pow(GS_norms[dim-1], 0.5)); //maximum possible value for x[dim-1]
	printf("max_num: %d\n", max_num);
	pthread_t threads[max_num];

	//create a lock for when each thread needs to edit shortest_vector
	pthread_mutex_t lock;
	if (pthread_mutex_init(&lock, NULL) != 0) {
		printf("Error: Mutex initialization failed\n");
		exit(1);
  }
	
	//divide the enumeration into threads by x[dim-1] value
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
