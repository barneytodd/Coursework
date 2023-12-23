#include <string.h>

double InnerProduct(int dim, int *arr1, int *arr2) {
  double sum1;
  int i;
  for (i=0, i<dim, i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sqrt(sum1);
}

bool EndCheck(int dim, int current_vec[dim], double (*A)[dim]) {
  //no minuses of ones already there
  //if all possible next vectors have dot product > 90 degrees, end
  //no repeats of what we've just tried
  //possibly have a dot product lookup_table
  bool possible_vec_check[dim][2];
  for (i=0, i<dim, i++) {
    if (current_vec[i] > 0) {
      possible_vec_check[i][0] = 1;
    else if (current_vec[i] < 0) {
      possible_vec_check[i][1] = 1;
    }
    }
  }
}

double ShortestVector(int dim, double (*A)[dim]) { //A is  amtrix of row vectors stacked
  //create an int array to keep track of current vector
  //define end points
  //do an initial branch
  //enumerate out to in
  //keep track of shortest vector
  int current_vec[dim];
  
}
