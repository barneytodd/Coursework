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
    int i, j, k;
    double sum1;
    bool skip[2] = {FALSE, FALSE};
    double B[dim-1, dim-1];
    if (dim == 1) {
        return A[0][0]
    }
    for (i=0, i<dim, i++) {
        skip[0] = skip[1] = FALSE;
        for (j=0; j<dim; j++) {
            if (j==i) {
                skip[0] = TRUE;
            }
            for (k=0; k<dim; k++) {
                if (k==i) {
                    skip[1] = TRUE;
                }
                B[j][k] = A[j+skip[0]][k+skip[1]]
            }
        }
        sum1 += pow(-1, i) * A[0][i] * Determinant(dim-1, B);
    }
}

double LimitCalc(int dim, double (*A)[dim]) {
    double gamma = tgamma(dim/2 + 1);
    double det = Determinant(A);
    return 1.05*(pow(gamma, 1/dim)/sqrt(M_PI))*pow(det, 1/dim);
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
        limit = LimitCalc(dim, A);
        printf("For Dimension: %d Limit: %.4f Got: %.4f\n", dim, limit, act);
        if(act > limit)
            printf("Limit %.4f, got %.4f\n", limit, act);
        assert(act <= expr);
    }
}

int main() {
    runTests(5, 80);
    runTests(7, 216);
    runTests(20, 114624);
    runTests(30, 14098308);
    return 0;
