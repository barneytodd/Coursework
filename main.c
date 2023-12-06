#include "LLL_Reduction.h"

int main() {
  double A[3][3];
  double B[3][3];
  LLL(3, 3, 0.75, A, B, (int[]){1, 1, 1}, (int[]){-1, 0, 2}, (int[]){3, 5, 6});
  printf("Orthonormalized Vectors (A):\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
  return 0;
}
