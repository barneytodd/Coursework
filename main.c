#include "LLL_Reduction.h"
#include "Enumeration.h"
//#include "test.h"
#include <stdio.h>
#include <stdlib.h>


#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdarg.h>


double SumArray(int dim, double arr[dim]) {
    int i;
    double sum1;
    for (i=0; i<dim; i++) {
        sum1 += arr[i];
    }
    return sum1;
}

double Determinant(int dim, double (*A)[dim]) {
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

double LimitCalc(int dim, double (*A)[dim]) {
    double gamma = tgamma(dim/2 + 1);
    printf("gamma %.4f\n", gamma);
    double det = Determinant(dim, A);
    printf("det: %.4f\n", det);
    return 1.05*(pow(gamma, 1/dim)/sqrt(M_PI))*pow(det, 1/dim);
}

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
  double A[dim][dim];
  for (i = 0; i < dim; i++) {
    for (j=0; j < dim; j++) {
      k = 1 + dim*i + j;
      //printf("argv[k]: %s\n", argv[k]);
      if (argv[k][0] == '[') {
        //printf("argv[k][1]: %.4f\n", strtod(&argv[k][1], NULL));
        A[i][j] = strtod(&argv[k][1], NULL); //Null might need to be changed here
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
  double limit = LimitCalc(dim, A);
  printf("Limit: %.4f\n", limit);
      
  //double vec1[3] = {1.0, 1.0, 1.0};
  //double vec2[3] = {-1.0, 0.0, 2.0};
  //double vec3[3] = {3.0, 5.0, 6.0};
  LLL(0.75, dim, A);
  printf("Orthonormalized Vectors (A):\n");
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
  //printf("Maximum value of double: %e\n", DBL_MAX);
  double shortest_length = ShortestVector(dim, A);
  FILE *result = fopen("result.txt", "w");
  fprintf(result, "%.4f\n", shortest_length);
  return 0;
}
