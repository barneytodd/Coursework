#include "LLL_Reduction.h"
#include "Enumeration.h"
#include "GeneralFunctions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv) {
  int dim, i, j, k;
	
	//calculate the input dimnesion by looking for the start of the second vector
  if (argc>2) {
    for (i = 2; i < argc; i++) {
      if (argv[i][0] == '[') {
        dim = i-1;
        break;
      }
    }
  }
  else if (argc == 2) {
    dim = 1;
  }
	else {
		printf("Error: no input vectors provided\n");
		exit(1);
	}
  
  //check that there are dim^2 + 1 arguments
  if (dim != (int)pow(argc-1, 0.5)) {
    printf("Error: Incorrect input format\nDimension of first vector: (%d) should equal sqrt(num input arguments): (%d)\n", dim, (int)pow(argc-1, 0.5));
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
      FreeMatrix(i, &A);
			perror("Failed to allocate memory for the rows of the input matrix");
			exit(1);
    }
  }

  //load the input vectors into A, and check for incorrect input formats
  char *endptr;
  for (i = 0; i < dim; i++) {
    for (j=0; j < dim; j++) {
      k = 1 + dim*i + j;
      //check that each vector starts with '[number'
      if (j == 0) {
        if (argv[k][0] != '[') {
          printf("Error: Incorrect input format\nExpected format for start of vector: '[number'\nInput format: '%s'\n", argv[1]);
					FreeMatrix(dim, &A);
          exit(1);
        }
        A[i][j] = strtod(&argv[k][1], &endptr); 
        //check the format of 1 dimensional inputs
        if (dim==1 && strcmp(endptr, "]") != 0) {
          printf("Error: Incorrect input format\nDimension = 1\nExpected format: '[number]'\nInput format: '%s'\n", argv[1]);
					FreeMatrix(dim, &A);	
          exit(1);
        }
        else if (dim==1) {
          continue;
        }
      }
      //check that the each vector has no more than dim elements and that the last element is in the format 'number]'
      else if (j == dim-1) {
        A[i][j] = strtod(argv[k], &endptr);
        if (strcmp(endptr, "]") == 0) {
          continue;
        }
        else if (strcmp(endptr, "") == 0) {
          printf("Error: Incorrect input format\nVector %d has too many elements\nExpected: %d elements\n", i+1, dim);
					FreeMatrix(dim, &A);	
          exit(1);
        }
        else {
          printf("Error: Incorrect input format\nExpected format for end of vector: 'number]'\nInput format: '%s'\n", argv[k]);
					FreeMatrix(dim, &A);		
          exit(1);
        }
      }
      //check that all the other numbers are formatted correctly
      else {
        A[i][j] = strtod(argv[k], &endptr);
      }
      if (strcmp(endptr, "]") == 0) {
        printf("Error: Incorrect input format\nExpected square matrix\nVector %d is of length %d, should be length %d\n", i+1, j+1, dim);
				FreeMatrix(dim, &A);	
        exit(1);
      }
      else if (strcmp(endptr, "") != 0) {
        printf("Error: Incorrect input format\nExpected format for all but the last entry of each vector: 'number' or '[number'\nInput format: '%s'\n", argv[k]);
				FreeMatrix(dim, &A);	
        exit(1);
      }
    }
  }
  endptr = NULL;

	//if dim = 1, the shortest vector is A[0][0]
	if (dim==1) {
		FILE *result = fopen("result.txt", "w");
	  if (result == NULL) {
	    perror("Error opening the result file");
	    exit(1);
	  }
	  fprintf(result, "%.4f\n", A[0][0]);
	  if (fclose(result) != 0) {
	    perror("Error closing the result file");
	    exit(1);
  	}
		FreeMatrix(dim, &A);
		return 0;
	}

  //initialise the matrix B, which stores GS orthogonalised values
	double **B = (double **)malloc(dim * sizeof(double *)); 
	if (B == NULL) {
			perror("Failed to allocate memory for the B matrix");
			FreeMatrix(dim, &A);
			exit(1);
	}
	for (i=0; i<dim; i++) {
		B[i] = (double *)malloc(dim * sizeof(double));
		if (B[i] == NULL) {
			FreeMatrix(dim, &A);
			FreeMatrix(i, &B);
			perror("Failed to allocate memory for the rows of B");
			exit(1);
		}
	}

	//initialise the array Mu, which stores the mu values during GS orthogonalisation
	//Mu can be stored as a matrix, however only the values below the diagonal are useful, so we store it as an 1D array to save memory
	//Mu[(i-1)*i+j] in 1D array form is the equivalent of Mu[i][j] in matrix form
	double *Mu = (double *)malloc((dim-1)*dim/2 * sizeof(double *)); //stores Mu values for GramSchmidt orthogonalisation
	if (Mu == NULL) {
		FreeMatrix(dim, &A);
		FreeMatrix(dim, &B);
		perror("Failed to allocate memory for Mu");
		exit(1);
	}

  //reduce the lattice basis using Lenstra–Lenstra–Lovász lattice reduction
  //LLL(0.75, dim, A, B, Mu);

  //compute lattice enumeration to find the shortest vector
  double shortest_length = 922.0;//ShortestVector(dim, A, B, Mu);
  printf("shortest length: %.4f\n", shortest_length);
  //free the memory allocated for A
  FreeMatrix(dim, &A);
	FreeMatrix(dim, &B);
	free(Mu);
	Mu = NULL;		

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
