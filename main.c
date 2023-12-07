#include "LLL_Reduction.h"
#include <stdio.h>

int main() {
  double A[3][3];
  double B[3][3];
  LLL(3, 0.75, 3, A, B, (int[]){1, 0, 0}, (int[]){0, 1, 0}, (int[]){0, 0, 1});
  printf("Orthonormalized Vectors (A):\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
  return 0;
}
