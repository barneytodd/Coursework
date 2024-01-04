#ifndef LLL_REDUCTION_H
#define LLL_REDUCTION_H

void LLL(double delta, int dim, double **A);

void GramSchmidt(int dim, int start, double B[][dim]);

void update_matrices(int dim, int start, double **A, double B[][dim]);

double InnerProduct(int dim, double *arr1, double *arr2);

#endif
