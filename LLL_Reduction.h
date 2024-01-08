#ifndef LLL_REDUCTION_H
#define LLL_REDUCTION_H

#include <stdbool.h>

void LLL(double delta, int dim, double **A, double **B, double *Mu);

void GramSchmidt(int dim, int start, double **B, double *Mu);

void update_matrices(int dim, int start, double **A, double **B, double *Mu);

bool CheckOrth(int dim, int start, double **B);

#endif
