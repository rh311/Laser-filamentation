#include <complex.h>

void slice(double complex *input_array, double complex **output_array, int start, int step, int size) {
  double complex *ptr = input_array + start;
  for (int i = 0; i < size; ++i) {
	output_array[i] = ptr;
	ptr += step;
  }
}

void copy_vector(double *vector_a, double *vector_b, int n) {
  for (int i = 0; i < n; ++i)
	vector_a[i] = vector_b[i];
}

void multiply_vector(double *vector_a, double *vector_b, double number_c, int n) {
  if (number_c == 1.)
	copy_vector(vector_a, vector_b, n);
  else
	for (int i = 0; i < n; ++i)
	  vector_a[i] = number_c * vector_b[i];
}

