#include "LLL_Reduction.h"
#include <stdio.h>

int main() {
  double A[3][3];
  //double v1[3] = 
  //double v2[3] = 
  //double v3[3] = 
  double vec1[3] = {1.0, 1.0, 1.0};
  double vec2[3] = {-1.0, 0.0, 2.0};
  double vec3[3] = {3.0, 5.0, 6.0};
  LLL(0.75, 3, A, vec1, vec2, vec3);
  printf("Orthonormalized Vectors (A):\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
  return 0;
}
