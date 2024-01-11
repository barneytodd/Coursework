#include "LLL_Reduction.h"
#include "Enumeration.h"
#include "GeneralFunctions.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <time.h>


//used to determine whether the input matrix is an identity matrix, returns the sum of the non-diagonal entries of a row
double SumArray(int dim, int i, double *arr) { 
    int j;
    double sum1 = 0;
    for (j=0; j<i; j++) {
        sum1 += arr[j];
    }
    for (j=i+1; j<dim; j++) {
        sum1 += arr[j];
    }
    return sum1;
}


//returns the determinant of a matrix, using LU Decomposition
double Determinant(int dim, double **A, bool *check) { 
    int i, j, k;
    double U[dim][dim];
    double factor;
    double determinant = 1.0;
    for (i=0; i<dim; i++) {
        for (j=0;j<dim;j++) {
            U[i][j] = A[i][j];
        }
    }
    //calculates the upper triangular matrix U, where A = L x U 
		//(L is a lower triangular matrix with 1's on the leading diagonal, therfore its determinant is 1 and so det(A) = det(U))
    for (i=0;i<dim;i++) {
        for (j=i+1;j<dim;j++) {
            factor = U[j][i]/U[i][i];
            for (k=i;k<dim;k++) {
                U[j][k] -= U[i][k]*factor;
            }
        }
        //if determinant gets too big, we can compute a running calculation of det^(1/dim) instead
        //best not to do this unless we have to, as it increases the inaccuracy caused by inaccuracy in calculations with doubles
        if (*check) {
            determinant *= pow(fabs(U[i][i]), 1.0/dim);
        }
        else {
            if (isinf(determinant*U[i][i])) {
                determinant = pow(fabs(determinant), 1.0/dim);
                determinant *= pow(fabs(U[i][i]), 1.0/dim);
                *check = true;
            }
            else {
                determinant *= U[i][i];
            }
        }
    }
    return determinant;
}


//returns an estimate for an upper bound of the shortest vector length, using the equation 1.05 * (gamma(n/2+1))^(1/dim)/sqrt(pi) * det(input matrix)^(1/dim)
double LimitCalc(int dim, double **A) { 
    bool det_check = false; //checks whether determinant has alredy been raised to the power of 1/dim
    double gamma = tgamma((float)dim/2 + 1);
    double det = Determinant(dim, A, &det_check);
    if (!det_check) {
        det = pow(fabs(det), 1.0/dim);
    }
    return 1.2*(pow(gamma, 1.0/dim)/sqrt(M_PI))*det;
}


//runs a test to check whether the shortest vector returned is a plausible solution or not
//Deals with two cases: (a) the input matrix is an identity matrix or (b) it is not
void runTests(int dim, double **A, double **B, double *Mu) { 
    
    int i;
    double limit = LimitCalc(dim, A);

    LLL(0.75, dim, A, B, Mu);
    double shortest_vector = ShortestVector(dim, A, B, Mu);

    //true if the input matrix is an identity matrix
    bool unit_test = true;
    for (i=0; i<dim; i++) {
        if (A[i][i] != 1.0 || SumArray(dim, i, A[i]) != 0.0) { 
            unit_test = false;
            break;
        }
    }

    //if the input matrix is an identity matrix, we expect to get the shortest vector length to be 1.0
    if (unit_test) {
        printf("For Dimension: %d Expected: %.4f Got: %.4f\n", dim, 1.0, shortest_vector);
        if(shortest_vector != 1.0) {
            FreeMatrix(dim, &A);
						FreeMatrix(dim, &B);
						free(Mu);
						Mu = NULL;
        }
        assert(shortest_vector == 1.0);
    }
        

    //if the input matrix is not an identity matrix, we expect to get a shortest vector length below the estimated upper bound   
    else {
        printf("For Dimension: %d Upper bound estimate: %.4f Got: %.4f\n", dim, limit, shortest_vector);
        if(shortest_vector > limit) {
					printf("Limit only an estimate of the upper bound, if shortest_vector is close to the limit, it may still be the correct answer\n");
					FreeMatrix(dim, &A);
					FreeMatrix(dim, &B);
					free(Mu);
					Mu = NULL;
        }
        assert(shortest_vector <= limit);
    }
}

//calls runTests, first with an identity matrix, and then with a randomly generated matrix
int main() {
	int i, j;
	int dim = 100; //dimension of first matrix

	//initialise the input matrix A
	double **A = (double **)calloc(dim, sizeof(double *));
	if (A==NULL) {
			perror("failed to allocate memory for the input matrix");
			exit(1);
	}
	for (i=0; i<dim; i++) {
			A[i] = (double *)calloc(dim, sizeof(double));
			if (A[i]==NULL) {
					FreeMatrix(i, &A);
					perror("failed to allocate memory for the rows of the input matrix");
					exit(1);
			}                
	}
	double **B = (double **)malloc(dim * sizeof(double *)); //stores GS orthogonalised values
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
				perror("Failed to allocate memory for the rows of the input matrix");
				exit(1);
		}
	}
    
	double *Mu = (double *)malloc((dim-1)*dim/2 * sizeof(double)); //stores Mu values for GramSchmidt orthogonalisation
	if (Mu == NULL) {
	FreeMatrix(dim, &A);
	FreeMatrix(dim, &B);
	perror("Failed to allocate memory for Mu");
	exit(1);
	}
    
  
	//set A to be the dim x dim identity matrix
	for (i=0; i<dim; i++) { 
			A[i][i] = 1.0;
	}
	
	printf("Input matrix dim: %d, identity matrix\n", dim);
	runTests(dim, A, B, Mu); 

	//set the bounds for the values of the second input matrix, and its dimension
	int min = -10000;
	int max = 10000;
	int new_dim = 10;
	
	//free any rows no longer needed
	if (new_dim<dim) {
		for (i=new_dim; i<dim; i++)  {
			free(A[i]);
			free(B[i]);
			A[i] = B[i] = NULL;
		}
	}

	dim = new_dim;
	//resize A, B and Mu
	A = realloc(A, dim * sizeof(double *));
	if (A==NULL) {
			perror("Failed to reallocate memory for the input matrix");
			FreeMatrix(dim, &B);
			free(Mu);
			Mu = NULL;
			exit(1);
	}			
	for (i=0;i<dim;i++) {
			A[i] = realloc(A[i], dim * sizeof(double));
			if (A[i]==NULL) {
					FreeMatrix(i, &A);
					FreeMatrix(dim, &B);
					free(Mu);
					Mu = NULL;
					perror("Failed to reallocate memory for the rows of the input matrix");
					exit(1);
			}          
	}
    
	B = realloc(B, dim * sizeof(double *)); //stores GS orthogonalised values
	if (B == NULL) {
			perror("Failed to reallocate memory for the B matrix");
			FreeMatrix(dim, &A);
			free(Mu);
			Mu = NULL;
			exit(1);
	}
	for (i=0; i<dim; i++) {
		B[i] = realloc(B[i], dim * sizeof(double));
		if (B[i] == NULL) {
				FreeMatrix(dim, &A);
				FreeMatrix(i, &B);
				free(Mu);
				Mu = NULL;
				perror("Failed to reallocate memory for the rows of the input matrix");
				exit(1);
		}
	}
		
	Mu = realloc(Mu, (dim-1)*dim/2 * sizeof(double)); //stores Mu values for GramSchmidt orthogonalisation
	if (Mu == NULL) {
		FreeMatrix(dim, &A);
		FreeMatrix(dim, &B);
		perror("Failed to reallocate memory for Mu");
		exit(1);
	}
	
	//initialise A to a set of random doubles, sampled from U(min, max)
	srand(time(NULL));
	for (i=0;i<dim;i++) {
			for (j=0;j<dim;j++) {
					A[i][j] = ((double)rand() / RAND_MAX) * (max-min) + min; 
			}
	}

	printf("Input matrix dim: %d, randomly generated matrix with values between %d and %d\n", dim, min, max);
	runTests(dim, A, B, Mu);

	FreeMatrix(dim, &A);
	FreeMatrix(dim, &B);
	free(Mu);
	Mu = NULL;
	return 0;
}
