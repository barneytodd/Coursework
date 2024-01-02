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

double Determinant(int dim, double **A) {
    int i, j, k;
    double sum1;
    bool skip[2] = {false, false};
    double B[dim-1][dim-1];
    if (dim == 1) {
        return A[0][0];
    }
    for (i=0; i<dim; i++) {
        skip[0] = skip[1] = false;
        for (j=0; j<dim; j++) {
            if (j==i) {
                skip[0] = true;
            }
            for (k=0; k<dim; k++) {
                if (k==i) {
                    skip[1] = true;
                }
                B[j][k] = A[j+skip[0]][k+skip[1]];
            }
        }
        sum1 += pow(-1, i) * A[0][i] * Determinant(dim-1, B);
    }
    return sum1;
}

double LimitCalc(int dim, double **A) {
    double gamma = tgamma(dim/2 + 1);
    double det = Determinant(dim, A);
    return 1.05*(pow(gamma, 1/dim)/sqrt(M_PI))*pow(det, 1/dim);
}

void runTests(int dim, double **A)
{
    //typedef struct {
    //    double *array;
    //    size_t size;
    //} Array;
    
    int i, j;
    //double A[dim][dim];
    //double *arr;
    //double limit;
    
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

    //else {
    //    limit = LimitCalc(dim, A);
    //    printf("For Dimension: %d Limit: %.4f Got: %.4f\n", dim, limit, shortest_vector);
    //    if(shortest_vector > limit) {
    //        printf("Limit %.4f, got %.4f\n", limit, shortest_vector);
    //    }
    //    assert(shortest_vector <= limit);
    //}
    //va_end(args);
}

int main() {
    int i;
    int dim = 5;
    double **A = (double **)calloc(dim, sizeof(double *));
    if (A==NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }
    
    for (i=0; i<dim; i++) {
        A[i] = (double *)calloc(dim, sizeof(double));
        if (A[i]==NULL) {
            printf("Memory allocation failed\n");
            return 1;
        }                
    }

    for (i=0; i<dim; i++) {
        A[i][i] = 1.0;
    }
    
    runTests(dim, A); 
    //runTests(7, 216);
    //runTests(20, 114624);
    //runTests(30, 14098308);
    return 0;
}
