void FreeMemoryA(int dim, double **A) {
	int i;
	for (i=0; i<dim; i++) {
		free(A[i]);
		A[i] = NULL;
	}
	free(A);
	A = NULL;
}

void FreeMemoryMu(int dim, double *Mu) {
	free(Mu);
	B = NULL;
}

void FreeMemory(int dim, int *choose_matrix, ...) { //double **A, double **B, double *Mu) {
	va_list ap;
	va_start (ap, choose_matrix);
	while (*choose_matrix > 0) {
		switch(*choose_matrix) {
			case 1:
				FreeMemoryA(dim, va_arg(ap, double **));
				choose_matrix++;
				break;
			case 2:
				FreeMemoryMu(dim, va_arg(ap, double *));
				choose_matrix++;
				break;
		}
	}
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
