#include <math.h>
#include "math_functions.h"
#include "physical_expressions.h"

extern double pi;
extern double speedOfLight;

double n0(double wavelength) {
  return sqrt(1. + 0.92549 / (1. - square(0.000007376 / wavelength)) + 6.96747 / (1. - square(0.003279 / wavelength)));
}

double GVD_func(double wavelength) {
  return (cube(wavelength) / (2. * pi * square(speedOfLight))) * second_derivative(n0,
																				   wavelength);
}

double z_sf(double omega, double n0, double rad_r, double PPcr, double focal_length) {
  double z_sf_0 = (0.367 * omega * n0 * square(rad_r) / speedOfLight) / sqrt(square(sqrt(PPcr) - 0.852) - 0.0217);
  return z_sf_0 * focal_length / (z_sf_0 + focal_length);
}

double Power_cr(double wavelength, double n0, double n2_I) {
  return 3.77 * square(wavelength) / (8. * pi * n0 * n2_I);
}


