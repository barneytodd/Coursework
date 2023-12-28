#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

double SumArray(int dim, double arr[dim]) {
    int i;
    double sum1
    for (i=0; i<dim; i++) {
        sum1 += arr[i];
    }
    return sum1
}

double Determinant(int dim, double (*A)[dim]) {
    
}

double LimitCalc(int dim, double (*A)[dim]) {
    double gamma = tgamma(dim/2 + 1);
    double det;
    
}

void runTests(int dim, ...)
{
    //typedef struct {
    //    double *array;
    //    size_t size;
    //} Array;
    
    int i, j, k;
    double A[dim][dim];
    double *arr;
    double limit;
    
    va_list args;
    va_start (args, dim);
    for (i = 0; i < dim; i++) {
        arr = va_arg(args, double*);
        for (j=0; j < dim; j++) {
            A[i][j] = arr[j];
        }
    }
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
    double shortest_length = ShortestVector(3, A);

    bool unit_test = TRUE;
    for (i=0, i<dim, i++) {
        if (SumArray(dim, A[i]) != 1.0) {
            unit_test = FALSE;
            break
        }
    }

    if (unit_test) {
        printf("For numbers: %d Expected: %llu Got: %llu\n", n, expr, act);
        if(shortest_vector != 1.0) {
            printf("Expected %llu, got %llu\n", expr, act);
        }
        assert(act == 1.0);
    }

    else {
        limit = 
        printf("For numbers: %d Expected: %llu Got: %llu\n", n, expr, act);
        if(act > expr)
            printf("Expected %llu, got %llu\n", expr, act);
        assert(act <= expr);
    }
}

int main() {
    runTests(5, 80);
    runTests(7, 216);
    runTests(20, 114624);
    runTests(30, 14098308);
    return 0;
