#include <string.h>

double InnerProduct(int dim, int *arr1, int *arr2) {
  double sum1;
  int i;
  for (i=0, i<dim, i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sqrt(sum1);
}

bool DotProductCheck(int dim, int current_vec[dim], bool possible_vec_check[dim][2], double (*A)[dim]) {
  //no minuses of ones already there
  //if all possible next vectors have dot product > 90 degrees, end
  //no repeats of what we've just tried
  //possibly have a dot product lookup_table

  int i, j;
  double sum1;
  double norm = InnerProduct(dim, current_vec, current_vec);;
  //check dot products
  for (i=0; i<dim; i++) {
    sum1 = 0;
    norm *= InnerProduct(dim, A[i], A[i]);
    for (j=0; j<dim; j++) {
      sum1 += InnerProduct(dim, A[i], current_vec[j]*A[j]);
    }
    sum1 /= norm;
    if (sum1 > 0) {
      return True;
    }
  }
  return False;
}

double ShortestVector(int dim, double (*A)[dim]) { //A is  amtrix of row vectors stacked
  //create an int array to keep track of current vector
  //define end points
  //do an initial branch
  //enumerate out to in
  //keep track of shortest vector
  int current_vec[dim];
  int last_vec;

  //check possible vecs using rule 1
  bool possible_vec_check[dim][2];
  for (i=0; i<dim; i++) {
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
