#include "LLL_Reduction.h"
#include "Enumeration.h"
#include <stdio.h>
#include <stdlib.h>

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
  //double shortest_length = ShortestVector(3, A);
  //FILE *result = fopen("result.txt", "w");
  //fprintf(result, "%.4f\n", shortest_length);
  return 0;
}
