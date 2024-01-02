#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdarg.h>
#include "LLL_Reduction.h"
#include "Enumeration.h"
#include <stdbool.h>


double SumArray(int dim, double arr[dim]) {
    int i;
    double sum1;
    for (i=0; i<dim; i++) {
        sum1 += arr[i];
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
    for (i=0;i<dim-1;i++) {
        for (j=i+1;j<dim;j++) {
            factor = A[i][j]/A[i][i];
            for (k=i;k<dim;k++) {
                U[j][k] -= U[i][k]*factor;
            }
        }
        determinant *= U[i][i];
    }
    printf("determinant: %.4f\n", determinant);
    return determinant;
}

double LimitCalc(int dim, double **A) {
    double gamma = tgamma(dim/2 + 1);
    printf("gamma: %.4f\n", gamma);
    double det = Determinant(dim, A);
    return 1.05*(pow(gamma, 1/dim)/sqrt(M_PI))*pow(det, 1/dim);
}

void runTests(int dim, double **A)
{
    
    
    int i, j;
    //double A[dim][dim];
    //double *arr;
    double limit;
    
    //va_list args;
    //va_start (args, dim);
    //for (i = 0; i < dim; i++) {
    //    arr = va_arg(args, double*);
    //    for (j=0; j < dim; j++) {
    //        A[i][j] = arr[j];
    //    }
    //}
    
    printf("A\n");
    for (i = 0; i < dim; i++) {
        for (j=0; j < dim; j++) {
          printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
    
    LLL(0.75, 3, A);
    printf("Orthonormalized Vectors (A):\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
    double shortest_vector = ShortestVector(3, A);

    bool unit_test = true;
    for (i=0; i<dim; i++) {
        if (SumArray(dim, A[i]) != 1.0) { //needs to be more specific
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
    int dim = 5;
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
    //runTests(7, 216);
    //runTests(20, 114624);
    //runTests(30, 14098308);
    for (i=0;i<dim;i++) {
        free(A[i]);
    }
    free(A);
    return 0;
}
