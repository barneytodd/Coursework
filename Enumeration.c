double InnerProduct(int dim, int *arr1, int *arr2) {
  double sum1;
  int i;
  for (i=0; i<dim; i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sqrt(sum1);
}

double ShortestVector(int dim, double (*A)[dim]) {
	int i, j;
	double shortest_vector = InnerProduct(dim, A[0], A[0]);
	double current_norm;
	for (i=1; i<dim; i++) {
		current_norm = InnerProduct(dim, A[i], A[i]);
		if (current_norm < shortest_vector) {
			shortest_vector = current_norm;
		}
	}
	double Mu[dim][dim];
	double GS_norms[dim]; //may not need, as we can just change A
	for (i=1; i<dim; i++) {
		GS_norms[i-1] = InnerProduct(dim, A[i-1], A[i-1]);
		for (j=0; j<i; j++) {
			Mu[i][j] = InnerProduct(dim, A[i], A[j])/GS_norms[j];
		}
		A[i] = 
	}
}
