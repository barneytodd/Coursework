#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdarg.h>
#include "LLL_Reduction.h"
#include "Enumeration.h"
#include "Threads.h"
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
double Determinant(int dim, double **A) { 
    int i, j, k;
    double U[dim][dim];
    double factor;
    double determinant = 1.0;
    for (i=0; i<dim; i++) {
        for (j=0;j<dim;j++) {
            U[i][j] = A[i][j];
        }
    }
    
    for (i=0;i<dim;i++) {
        for (j=i+1;j<dim;j++) {
            factor = U[j][i]/U[i][i];
            for (k=i;k<dim;k++) {
                U[j][k] -= U[i][k]*factor;
            }
        }
        determinant *= pow(fabs(U[i][i]), 1.0/dim);
    }
    printf("determinant: %.4f\n", determinant);
    return determinant;
}


//returns an estimate for an upper bound of the shortest vector length, using the equation 1.05 * (gamma(n/2+1))^(1/dim)/sqrt(pi) * det(input matrix)^(1/dim)
double LimitCalc(int dim, double **A) { 
    double gamma = tgamma((float)dim/2 + 1);
    double det = fabs(Determinant(dim, A));
    return 1.05*(pow(gamma, 1.0/dim)/sqrt(M_PI))*pow(det, 1.0/dim);
}


//runs a test to check whether the shortest vector returned is a plausible solution or not
//Deals with two cases: (a) the input matrix is an identity matrix or (b) it is not
void runTests(int dim, double **A) { 
    
    int i;
    double limit = LimitCalc(dim, A);
    
    LLL(0.75, dim, A);
    double shortest_vector = ShortestVector1(dim, A);

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
            printf("Expected %.4f, got %.4f\n", 1.0, shortest_vector);
        }
        assert(shortest_vector == 1.0);
    }
        

    //if the input matrix is not an identity matrix, we expect to get a shortest vector length below the estimated upper bound   
    else {
        printf("For Dimension: %d Upper bound estimate: %.4f Got: %.4f\n", dim, limit, shortest_vector);
        if(shortest_vector > limit) {
            printf("Upper bound estimate %.4f, got %.4f\n", limit, shortest_vector);
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
            for (j=0;j<dim;j++) {
                free(A[j]);
            }
            free(A);
            perror("failed to allocate memory for the rows of the input matrix");
            exit(1);
        }                
    }
    
    //set A to be the dim x dim identity matrix
    for (i=0; i<dim; i++) { 
        A[i][i] = 1.0;
    }
    
    printf("Input matrix dim: %d, identity matrix\n", dim);
    runTests(dim, A); 

    //set the bounds for the values of the second matrix, and its dimension
    int min = -10000;
    int max = 10000;
    dim = 80;

    //resize A
    A = realloc(A, dim * sizeof(double *));
    if (A==NULL) {
        perror("Failed to reallocate memory for the input matrix");
        exit(1);
    }
    for (i=0;i<dim;i++) {
        A[i] = (double *)realloc(A[i], dim * sizeof(double));
        if (A[i]==NULL) {
            for (j=0;j<i;j++) {
                free(A[j]);
            }
            free(A);
            perror("Failed to reallocate memory for the rows of the input matrix");
            exit(1);
        }          
    }
    srand(time(NULL));
    
    //initialise A to a set of random doubles, sampled from U(min, max)
    for (i=0;i<dim;i++) {
        for (j=0;j<dim;j++) {
            A[i][j] = ((double)rand() / RAND_MAX) * (max-min) + min; 
        }
    }

    printf("Input matrix dim: %d, randomly generated matrix with values between %d and %d\n", dim, min, max);
    runTests(dim, A);

    
    for (i=0;i<dim;i++) {
        free(A[i]);
        A[i] = NULL;
    }
    free(A);
    A = NULL;
    return 0;
}
