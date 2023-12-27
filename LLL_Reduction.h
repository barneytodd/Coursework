void LLL(double delta, int dim, double (*A)[dim]);

void GramSchmidt(int dim, int start, double B[][dim]);

void update_matrices(int dim, int start, double (*A)[dim], double B[][dim]);

double InnerProduct(int dim, double *arr1, double *arr2);
