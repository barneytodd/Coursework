double VectorNorm(int dim, int *arr) {
  double sum1;
  int i;
  for (i=0, i<dim, i++) {
    sum1 += arr[i]*arr[i];
  }
  return sqrt(sum1);
}

double ShortestVector(int dim, double (*A)[dim]) {
  bool boolarray[dim]; //checks if each vector has been changed
  for (i=0, i<dim, i++) {
    boolarray[i] = True; //initialises boolarray to all True
  }
  double temp_vector[dim];
  for (i=0, i<dim, i++) {
    for (j=i, j<dim, j++) {
      if (i!=j && boolarray[j]) {
        double i_length = VectorNorm(dim, A[][i]);
        double j_lenght;
        
      }
    }
  }
}
  
