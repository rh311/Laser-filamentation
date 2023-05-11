#ifndef FILAMENTATION_REFACTORING_PHYSICAL_FORMULAE_H
#define FILAMENTATION_REFACTORING_PHYSICAL_FORMULAE_H

double n0(double wavelength);
double GVD_func(double wavelength);
double z_sf(double omega, double n0, double rad_r, double PPcr, double focal_length);
double Power_cr(double wavelength, double n0, double n2_I);

#endif
