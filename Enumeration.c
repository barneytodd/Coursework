//edit DotProduct to  only include valid vectors from possible_bool_check

#include <string.h>

double InnerProduct(int dim, int *arr1, int *arr2) {
  double sum1;
  int i;
  for (i=0, i<dim, i++) {
    sum1 += arr1[i]*arr2[i];
  }
  return sqrt(sum1);
}

bool DotProductCheck(int dim, int current_vec_constits[dim], double current_vec[dim], bool possible_vec_check[dim][2], double (*A)[dim]) {
  //no minuses of ones already there
  //if all possible next vectors have dot product > 90 degrees, end
  //no repeats of what we've just tried
  //possibly have a dot product lookup_table

  int i, j, k;
  double sum1;
  double norm = InnerProduct(dim, current_vec, current_vec);;
  //check dot products
  for (i=0; i<dim; i++) {
    //for (k=0, k<2, k++)
    sum1 = 0;
    norm *= InnerProduct(dim, A[i], A[i]);
    for (j=0; j<dim; j++) {
      //if (possible_vec_check[j][0]) {
      sum1 += InnerProduct(dim, A[i], current_vec_constits[j]*A[j]);
      //}
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
  int current_vec_constits[dim]; //array of the constituant vectors that make up current_vec
  memset(current_vec_constits, 0, dim*sizeof(int)); //initialise to all 0s
  current_vec_constits[0] = 1; //start with the first vector
  double current_vec[dim]; //actual current vec
  memcpy(current_vec, A[0], dim*sizeof(double));

  //check possible vecs using rule 1
  bool possible_vec_check[dim][2]; //e.g. possible_vec_check[0][0] = True means the positive first vector is allowed
  for (i=0; i<dim; i++) {
    //initialise bool array to True
    possible_vec_check[i][0] = True;
    possible_vec_check[i][1] = True;

  possible_vec_check[0][0] = False

    //no minuses of ones alredy there
  //  if (current_vec_constits[i] > 0) {
    //  possible_vec_check[i][0] = False;
    //}
    //else if (current_vec_constits[i] < 0) {
     // possible_vec_check[i][1] = False;
    //} 
  }

  bool failed_vectors[dim];
	bool end_branch;
	int starting_vector = 0;
	bool end_branch = False;
		while (!end_branch) {
			while (DotProductCheck(dim, current_vec_constits, current_vec, possible_vec_check, A)) {
				for (i=starting_vector; i<dim; i++) {
					if (possible_vectors[i][1]) && (!failed_vectors[i]) {
						current_vector += A[i];
						current_vector_constit[i] += 1;
						possible_vectors[i][0] = False;
						break;
					}
					if (possible_vectors[i][0]) && (!failed_vectors[i]) {
						current_vector -= A[i];
						current_vector_constit[i] -= 1;
						possible_vectors[i][1] = False;
						break;
					}
					
				}
			}
			for (i=0; i<dim; i++) {
				failed_vectors[i] = False;
			}
			for (i=dim-1; i>=0; i--) {
				if (i==starting_vector) {
					end_branch = True;
					starting_vector += 1;
					break;
				}
				if (current_vector_constits[i] < 0) {
					current_vector += A[i] //maybe need to iterate;
					current_vector_constits[i] += 1;
					failed_vectors[i] = True;
					break;
				}
				else if (current_vector_constits[i] > 0) {
					current_vector -= A[i] //maybe need to iterate;
					current_vector_constits[i] -= 1;
					failed_vectors[i] = True;
					break;
				}
				
			}
		
		}
		current_vector = A[]
		for (i=0, i<dim, i++) {
			
		}	 
  
  
}
