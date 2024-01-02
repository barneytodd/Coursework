#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdarg.h>
#include "LLL_Reduction.h"
#include "Enumeration.h"
#include <stdbool.h>


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


double Determinant(int dim, double **A) { //Using LU Decomposition
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
            factor = A[j][i]/A[i][i];
            for (k=i;k<dim;k++) {
                U[j][k] -= U[i][k]*factor;
            }
        }
        determinant *= U[i][i];
    }
    printf("U\n");
    for (i=0;i<dim;i++) {
        for (j=0;j<dim;j++) {
            printf("%.4f\t", U[i][j]);
        }
        printf("\n");
    }
    printf("determinant: %.4f\n", determinant);
    return determinant;
}

double LimitCalc(int dim, double **A) {
    double gamma = tgamma((float)dim/2 + 1);
    double det = Determinant(dim, A);
    return 1.05*(pow(gamma, 1.0/dim)/sqrt(M_PI))*pow(det, 1.0/dim);
}

void runTests(int dim, double **A)
{
    
    
    int i, j;
    
    double limit;
    
    
    
    printf("A\n");
    for (i = 0; i < dim; i++) {
        for (j=0; j < dim; j++) {
          printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
    
    LLL(0.75, dim, A);
    printf("Orthonormalized Vectors (A):\n");
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
    double shortest_vector = ShortestVector(dim, A);

    bool unit_test = true;
    for (i=0; i<dim; i++) {
        //printf("A[i][i], SumArray: %.4f %.4f\n", A[i][i], SumArray(dim, i, A[i]));
        if (A[i][i] != 1.0 || SumArray(dim, i, A[i]) != 0.0) { //needs to be more specific
            unit_test = false;
            break;
        }
    }

    if (unit_test) {
        printf("For Dimension: %d Expected: %.4f Got: %.4f\n", dim, 1.0, shortest_vector);
        if(shortest_vector != 1.0) {
            printf("Expected %.4f, got %.4f\n", 1.0, shortest_vector);
        }
        assert(shortest_vector == 1.0);
    }

    else {
        limit = LimitCalc(dim, A);
        printf("For Dimension: %d Limit: %.4f Got: %.4f\n", dim, limit, shortest_vector);
        if(shortest_vector > limit) {
            printf("Limit %.4f, got %.4f\n", limit, shortest_vector);
        }
        assert(shortest_vector <= limit);
    }
    //va_end(args);
}

int main() {
    int i, j;
    int dim = 40;
    double **A = (double **)calloc(dim, sizeof(double *));
    if (A==NULL) {
        exit(1);
    }
    
    for (i=0; i<dim; i++) {
        A[i] = (double *)calloc(dim, sizeof(double));
        if (A[i]==NULL) {
            for (j=0;j<dim;j++) {
                free(A[j]);
            }
            free(A);
            exit(1);
        }                
    }

    for (i=0; i<dim; i++) {
        A[i][i] = 1.0;
    }
    printf("dim1: %d\n", dim);
    printf("A1: \n");
    for (i=0;i<dim;i++) {
        for (j=0;j<dim;j++) {
            printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
    runTests(dim, A); 

    int min = -10000;
    int max = 10000;
    dim = 10;
    A = realloc(A, dim * sizeof(double));

    
    
    for (i=0;i<dim;i++) {
        for (j=0;j<dim;j++) {
            A[i][j] = ((double)rand() / RAND_MAX) * (max-min) + min; //initialises A to random doubles sampled from Unif(min, max)
        }
    }
    printf("A2: \n");
    printf("[");
    for (i=0;i<dim;i++) {
        printf("[");
        for (j=0;j<dim;j++) {
            printf("%.4f, ", A[i][j]);
        }
        printf("] \n");
    }
    printf("] \n");
    runTests(dim, A);
    //runTests(20, 114624);
    //runTests(30, 14098308);
    for (i=0;i<dim;i++) {
        free(A[i]);
    }
    free(A);
    return 0;
}
