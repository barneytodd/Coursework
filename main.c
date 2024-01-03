#include "LLL_Reduction.h"
#include "Enumeration.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





int main(int argc, char *argv[]) {
  int dim, i, j, k;
  
  printf("You have entered %d arguments:\n", argc);
  for (i = 2; i < argc; i++) {
    if (argv[i][0] == '[' && i!=1) {
      dim = i-1;
      break;
    }
  }
  printf("dim: %d", dim);

  //initialise the input matrix A
  double **A = (double **)malloc(dim * sizeof(double *));
  if (A == NULL) {
      perror("Failed to allocate memory for the input matrix");
      exit(1);
  }
  for (i=0; i<dim; i++) {
    A[i] = (double *)malloc(dim * sizeof(double));
    if (A[i] == NULL) {
        for (j=0; j<i; j++) {
            free(A[j]);
        }
        free(A);
        perror("Failed to allocate memory for the rows of the input matrix");
        exit(1);
    }
  }

  //load the input vectors into A
  for (i = 0; i < dim; i++) {
    for (j=0; j < dim; j++) {
      k = 1 + dim*i + j;
      if (argv[k][0] == '[') {
        A[i][j] = strtod(&argv[k][1], NULL); 
        continue;
      }
      if (argv[k][-1] == ']') {
        A[i][j] = strtod(&argv[k][-1], NULL);
      }
      else {
        A[i][j] = strtod(argv[k], NULL);
      }
    }
  }
  
  printf("A\n");
  for (i = 0; i < dim; i++) {
    for (j=0; j < dim; j++) {
      printf("%.4f\t", A[i][j]);
    }
    printf("\n");
  }

  //reduce the lattice basis using Lenstra–Lenstra–Lovász lattice reduction
  LLL(0.75, dim, A);
  
  printf("Orthonormalized Vectors (A):\n");
  for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
          printf("%.4f\t", A[i][j]);
      }
      printf("\n");
  }

  //compute lattice enumeration to find the shortest vector
  double shortest_length = ShortestVector(dim, A);
  
  //free the memory allocated for A
  for (i=0;i<dim;i++) {
    free(A[i]);
  }
  free(A);

  //save the output to result.txt
  FILE *result = fopen("result.txt", "w");
  if (result == NULL) {
    perror("Error opening the result file");
    exit(1);
  }
  fprintf(result, "%.4f\n", shortest_length);
  if (fclose(result) != 0) {
    perror("Error closing the result file");
    exit(1);
  }
  
  return 0;
}
