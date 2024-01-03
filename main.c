#include "LLL_Reduction.h"
#include "Enumeration.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(int argc, char *argv[]) {
  int dim, i, j, k;

  //check if num rows = num columns
  //check if each row has the same num columns
  //check if there are spaces before/after the []
  //check if arguments are valid numbers
  printf("You have entered %d arguments:\n", argc);
  for (i = 2; i < argc; i++) {
    if (argv[i][0] == '[') {
      dim = i-1;
      break;
    }
  }
  printf("dim: %d", dim);

  //there should be dim^2 + 1 arguments
  if (dim != (int)pow(argc-1, 0.5)) {
    printf("Error: Incorrect input format\n");
    exit(1);
  }

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
    //if (dim == 1) {
      //if (argv[0][0] != '[' || argv[0][-1] != ']') {
        //printf("Error: Incorrect input format\n");
      //}
      //else {
      //  char *temp = &argv + 1
      //  A[0][0] = strtod(&argv[0][1]
      //}
    //}
    for (j=0; j < dim; j++) {
      char *endptr;
      k = 1 + dim*i + j;
      if (argv[k][0] == '[' && j == 0) {
        A[i][j] = strtod(&argv[k][1], &endptr); 
      }
      else if (argv[k][-1] == ']' && j == dim-1) {
        A[i][j] = strtod(argv[k], &endptr);
        //if (endptr == ']') {
        //  continue;
        //}
        //else {
        //  printf("Error: Incorrect input format\n");
        //}
      }
      else {
        A[i][j] = strtod(argv[k], &endptr);
      }
      //if (endptr != ' ') {
      //  printf("Error: Incorrect input format\n");
      //}
      printf("endptr: %s\n", endptr);
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
