#include <string.h>

double InnerProduct(int dim, int *arr1, int *arr2) {
  double sum1;
  int i;
  for (i=0, i<dim, i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sqrt(sum1);
}

int EndCheck(int dim, int current_vec[dim], int last_vec, double (*A)[dim]) {
  //no minuses of ones already there
  //if all possible next vectors have dot product > 90 degrees, end
  //no repeats of what we've just tried
  //possibly have a dot product lookup_table

  int i, j;
  double sum1;
  double norm;
  //check dot products
  bool dp_check = False
  for (i=0; i<dim; i++) {
    sum1 = 0;
    norm = InnerProduct(dim, A[i], A[i]);
    for (j=0; j<dim; j++) {
      sum1 += InnerProduct(dim, A[i], current_vec[j]*A[j]);
      norm *= InnerProduct(dim, A[j], A[j]);
    }
    sum1 /= norm;
    if (sum1 > 0) {
      dp_check = True;
      break
    }
  }
  if (dp_check == False) {
    return -1
  }
  
  bool possible_vec_check[dim][2];
  for (i=0; i<dim; i++) {
    //don't repeat last step
    if (i == last_vec) {
      possible_vec_check[i][0] = False;
      possible_vec_check[i][1] = False;
      continue
    }

    //initialise bool array to True
    possible_vec_check[i][0] = True;
    possible_vec_check[i][1] = True;

    //no minuses of ones alredy there
    if (current_vec[i] > 0) {
      possible_vec_check[i][0] = False;
    }
    else if (current_vec[i] < 0) {
      possible_vec_check[i][1] = False;
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
