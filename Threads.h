#ifndef THREADS_H
#define THREADS_H

#include <pthread.h>

struct ThreadArgs {
  int num;
  int dim;
  double *GS_norms;
  double *Mu;
  double *shortest_vector;
  int *max_num;
	pthread_mutex_t *lock;
};

void *Enumerate(void *args);

double ShortestVector(int dim, double **A);

#endif 