#ifndef FILAMENTATION_REFACTORING_MATH_FUNCTIONS_H
#define FILAMENTATION_REFACTORING_MATH_FUNCTIONS_H

#include <complex.h>

double square(double x);
double cube(double x);
double gaussian(double x, double mean, double width);
double cnorm(double complex z);
double second_derivative(double (*f)(double), double x);

#endif
