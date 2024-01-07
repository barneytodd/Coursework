void FreeMemory(int dim, double **A, double **B, double *Mu) {
	int i;
	for (i=0; i<dim; i++) {
		free(A[i]);
		A[i] = NULL;
		free(B[i]);
		B[i] = NULL;
	}
	free(A);
	A = NULL;
	free(B);
	B = NULL;	
	free(Mu);
	Mu = NULL;
}

void FreeMemoryA(int dim, double **A) {
	int i;
	for (i=0; i<dim; i++) {
		free(A[i]);
		A[i] = NULL;
	}
	free(A);
	A = NULL;
}

void FreeMemoryB(int dim, double **A) {
	int i;
	for (i=0; i<dim; i++) {
		free(B[i]);
		B[i] = NULL;
	}
	free(B);
	B = NULL;
}


//compute the inner product between two vectors
double InnerProduct(int dim, double *arr1, double *arr2) {
  double sum1 = 0;
  int i;
  for (i=0; i<dim; i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sum1;
}
