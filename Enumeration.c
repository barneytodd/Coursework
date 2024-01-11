#include "GeneralFunctions.h"
#include "Enumeration.h"
#include <math.h>
#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

// runs the lattice enumeration loop for each thread
void *Enumerate(void *args) {	
  struct ThreadArgs *thread_args = (struct ThreadArgs *)args;  // restructure the arguments into the form of the struct above
  
  int x[thread_args->dim], i, j, k, n;  // x stores the number of each basis vector used to reach each lattice point
  double l[thread_args->dim];  // stores the total contribution squared, of the combination of basis vectors stored in x, in the direction of the ith GS vector
  
  for (i=0; i<thread_args->dim-1; i++) {
    x[i] = 0;
  }

  // set x[dim-1] to the thread number, and start with i = dim1
  x[thread_args->dim-1] = thread_args->num;
  i = thread_args->dim-1;
  
  double sum2;  // stores the sum of x[j] * Mu[j][i] for j>i
  double sum3;  // stores the sum (j>i) of the l[j] values
  double m = 0.0;  // counts the number of iterations of the while loop, needs to be double in case num iterations > max int
  double max_its = pow(2.0, pow(thread_args->dim-1, 2));  // an upper bound of the number of iterations required
  
  double short_vec = *(thread_args->shortest_vector);  // keep the value of the shortest_vector, we may need to adjust if it gets changed by another thread
  
  // max_num may get updated by the other threads
  // in the case that max_num falls below num, we can exit this thread
  while (*thread_args->max_num >= thread_args->num) { 
    
    // calculate the l[j] values from i upwards
    for (j=thread_args->dim-1; j >= i; j--) { 
    sum2 = 0;
      for (k = j+1; k < thread_args->dim; k++) {
        sum2 += x[k] * (*(thread_args->Mu))[(k-1)*k/2+j]; 
      }
      l[j] = (x[j] + sum2) * (x[j] + sum2) * (*(thread_args->GS_norms))[j]; 	
    }
    
    // sum the l[j] values for j>=i
    sum3 = 0;
    for (j=i; j<thread_args->dim; j++) {
      sum3 += l[j];
    }

    if (sum3 < (*(thread_args->shortest_vector))*(*(thread_args->shortest_vector))) {
      // if i=0 and sum3 < (current shortest vector length)^2, we have a new shortest vector
      if (i==0) {
        if (sum3 != 0) {
          pthread_mutex_lock(thread_args->lock);
          if (sum3 < (*(thread_args->shortest_vector))*(*(thread_args->shortest_vector))) {  //check again whilst inside the lock, to make sure shortest_vector hasn't changed
            *(thread_args->shortest_vector) = sqrt(sum3);
            printf("shortest_vector: %.4f\n", *(thread_args->shortest_vector));
            *(thread_args->max_num) = floor(*(thread_args->shortest_vector)/pow((*(thread_args->GS_norms))[thread_args->dim-1], 0.5));
          }
          pthread_mutex_unlock(thread_args->lock);
        }
        x[0]++;
      }

      // if i != 1, subtract 1 from i and then set x[i] to be the minimum integer such that l[i] < shortest_vector^2 - sum3
      // if there is no such integer, add the 1 back to i and then add 1 to x[i]
      else {
        i--;
        sum2 = 0;
        for (k=i+1; k<thread_args->dim; k++) {
          sum2 += x[k] * (*(thread_args->Mu))[(k-1)*k/2+i]; 
        }
        x[i] = round(- sum2);  //the integer which minimises l[i], if this doesn't work then no other integer will
        l[i] = ((double)x[i] + sum2) * ((double)x[i] + sum2) * (*(thread_args->GS_norms))[i]; 
        
        if (l[i] < (*(thread_args->shortest_vector)) * (*(thread_args->shortest_vector)) - sum3) {
          // subtract 1 from x[i] until l[i] is no longer < shortest_vector^2 - sum3
          // then add 1 to x[i] to make x[i] the minimum possible integer such that l[i] < shortest_vector^2 - sum3
          n=0;
          do {
            x[i]--;
            sum2 = 0;
            for (k=i+1; k<thread_args->dim; k++) {
              sum2 += x[k] * (*(thread_args->Mu))[(k-1)*k/2+i];
            }
            l[i] = (x[i] + sum2) * (x[i] + sum2) * (*(thread_args->GS_norms))[i]; 
            n++;
            if (n>100) {
              printf("Error: thread %d failed\n", thread_args->num);
              FreeMatrix(thread_args->dim, thread_args->A);
              FreeMatrix(thread_args->dim, thread_args->B);
              free(*(thread_args->Mu));
              *(thread_args->Mu) = NULL;		
              free(*(thread_args->GS_norms));
              *(thread_args)-> GS_norms = NULL;
              exit(1);
            }
          } while (l[i] < (*(thread_args->shortest_vector)) * (*(thread_args->shortest_vector)) - sum3);
          
          x[i]++; 
				}

        else {
          i+=1;
          x[i]++;
        }
      }
    }
			
    // if sum3 > shortest_vector^2, increase i by 1 and then increase x[i] by 1
    else {
      // if shortest_vector has been changed by another thread, we need to perform some checks
      if (short_vec != *(thread_args->shortest_vector)) {
        short_vec = *(thread_args->shortest_vector);
        
        l[thread_args->dim-2] = pow(x[thread_args->dim-2] + x[thread_args->dim-1] * (*(thread_args->Mu))[(thread_args->dim-2)*(thread_args->dim-1)/2+thread_args->dim-2], 2) * (*(thread_args->GS_norms))[thread_args->dim-2]; 
        
        // if l[dim-2] + l[dim-1] < shortest_vector^2, then we are fine to carry on from i = dim-2
        if (l[thread_args->dim-2]+l[thread_args->dim-1] < pow(*(thread_args->shortest_vector), 2)) {
          i = thread_args->dim-2;
          continue;
        }
          
        // if l[dim-2] + l[dim-1] > shortest_vector^2, then x[dim-2] is now out of range w.r.t. the new shortest_vector
        // this could potentially lead to this thread terminating before it has checked all possible x values
        else {
          // if l[dim-2] calculated with x[dim-2] < l[dim-2] calculated with x[dim-2]-1, then x[dim-2] is below the new accepted range
          // therefore we haven't yet checked the x[dim2] values in the new accepted range, so we reset i to dim-1 and carry on
          if (l[thread_args->dim-2] < pow((x[thread_args->dim-2]-1) + (x[thread_args->dim-1]-1) * (*(thread_args->Mu))[(k-1)*k/2+j], 2) * (*(thread_args->GS_norms))[thread_args->dim-2]) {
            i = thread_args->dim-1;
            continue;
          }
            
          // in the opposite case, x[dim-2] is above the new accepted range, and so we have already checked all possibilities in this new range
          // therefore we terminate this thread
          else {
            break;
          }
        }
      }
      
      // otherwise we proceed as normal
      i++;
      // if we're trying to increment x[dim-1], we've reached the end of this thread
      if (i==thread_args->dim-1) {
        break;
      }
      x[i]++;
    }
    m+=1;
    if (m > max_its) {
      printf("Error: Enumeration loop for thread %d failed\n", thread_args->num);
      FreeMatrix(thread_args->dim, thread_args->A);
      FreeMatrix(thread_args->dim, thread_args->B);
      free(*(thread_args->Mu));
      *(thread_args->Mu) = NULL;		
      free(*(thread_args->GS_norms));
      *(thread_args)-> GS_norms = NULL;
      exit(1);
    }
  }
  pthread_exit(NULL);
}

// Enumerate the lattice to find the shortest vector
double ShortestVector(int dim, double **A, double **B, double *Mu) {
  int i, j;
  
  // check the size of A and B
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
		
  double shortest_vector = sqrt(InnerProduct(dim, A[0], A[0]));  // keeps track of current shortest vector
  double current_norm;  // stores the norm of each basis vector
  
  // find the shortest basis vector and set shortest_vector to be equal to its norm
  for (i=1; i<dim; i++) {
    current_norm = sqrt(InnerProduct(dim, A[i], A[i]));  // calculates norm of each vector in A, then compares to shortest vector
    if (current_norm < shortest_vector) {
      shortest_vector = current_norm;
    }
  }
  
  printf("shortest basis vector: %.4f\n", shortest_vector);
  
  // stores the squared norm of each GramSchidt orthogonalised vector
  double *GS_norms = (double *)malloc(dim * sizeof(double)); 
  if (GS_norms == NULL) {
    FreeMatrix(dim, &A);
    FreeMatrix(dim, &B);
    free(Mu);
    Mu = NULL;		
    exit(1);
  }

  for (i=0;i<dim;i++) {
    GS_norms[i] = InnerProduct(dim, B[i], B[i]);
  }

  int max_num = floor(shortest_vector/pow(GS_norms[dim-1], 0.5));  // maximum possible value for x[dim-1]
  
  // create a lock for when each thread needs to edit shortest_vector
  pthread_mutex_t lock;
  if (pthread_mutex_init(&lock, NULL) != 0) {
    printf("Error: Mutex initialization failed\n");
    FreeMatrix(dim, &A);
    FreeMatrix(dim, &B);
    free(Mu);
    Mu = NULL;		
    free(GS_norms);
    GS_norms = NULL;
    exit(1);
  }
  int batch_size = fmin(10, dim/2);  // max number of threads allowed to be open at one time
  int n = (max_num+1)/batch_size;  // number of full batches needed
  int m;  // determines how many threads we create in each iteration
  
  // divide the enumeration into threads by x[dim-1] value, maximum number of threads allowed at one time = batch_size
  for (i=0; i<n+1; i++) {
    if (max_num+1-batch_size*i >= batch_size) {
      m = batch_size;
    }
    else {
      m = (max_num+1) % batch_size;
    }
    pthread_t threads[m];
    struct ThreadArgs args[m];
    for (j=0; j<m; j++) {
      args[j].num = j + batch_size*i;
      args[j].dim = dim;
      args[j].GS_norms = &GS_norms;
      args[j].Mu = &Mu;
      args[j].shortest_vector = &shortest_vector;
      args[j].max_num = &max_num;
      args[j].lock = &lock;
      args[j].A = &A;
      args[j].B = &B;
      if (pthread_create(&threads[j], NULL, &Enumerate, (void *)&args[j]) != 0) {
        printf("Error creating thread %d\n", j);
        FreeMatrix(dim, &A);
        FreeMatrix(dim, &B);
        free(Mu);
        Mu = NULL;		
        free(GS_norms);
        GS_norms = NULL;
        exit(1);
      }
      else {
        printf("thread %d created\n", i*batch_size+j);
      }
    }
    // requires all threads to finish before moving on to the next batch
    for (j=0; j<m; j++) {
      pthread_join(threads[j], NULL);
    }
    n = (max_num+1)/batch_size;  // update in case max_num has changed
  }
  pthread_mutex_destroy(&lock);
  
  free(GS_norms);
  GS_norms = NULL;
  
  return shortest_vector;
}
