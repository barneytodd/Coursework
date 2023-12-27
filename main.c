#include "LLL_Reduction.h"
#include "Enumeration.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
  double A[3][3];
  
  printf("You have entered %d arguments:\n", argc);
  for (int i = 0; i < argc; i++) {
      printf("%s\n", argv[i]);
  }
  
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
  double shortest_length = ShortestVector(3, A);
  printf("shortest length: %.4f\n", shortest_length);
  return 0;
}
