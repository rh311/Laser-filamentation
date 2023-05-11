#include <math.h>
#include <complex.h>
#include "math_functions.h"

double square(double x) {
  return x * x;
}

double cube(double x) {
  return x * x * x;
}

double gaussian(double x, double mean, double width) {
  return exp(-square((x - mean) / width) / 2.);
}

double cnorm(complex double z) {
  return square(creal(z)) + square(cimag(z));
}

double second_derivative(double (*f)(double), double x) {
  double dx = 0.001;
  double result = 0.;
  double intermediate_result;
  int counter = 0;
  do {
	intermediate_result = result;
	result = ((*f)(x + dx) - 2. * (*f)(x) + (*f)(x - dx)) / square(dx);
	dx /= 2.;
	++counter;
  } while (fabs(intermediate_result - result) > 0.01 * fabs(result));
  return result;
}