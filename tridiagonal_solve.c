#include <complex.h>
#include <stddef.h>
#include "grid.h"
#include "math_functions.h"
#include "create_array.h"
#include "tridiagonal_solve.h"

static double complex *Utmp;

void tridiagonal_solve_init(int max_matrix_size) {
  Utmp = create_array(max_matrix_size * sizeof(double complex));
}

void set_r_laplacian_coefficients(struct grid_t_r_z grid, double *a_r_lapl, double *b_r_lapl, double *c_r_lapl) {
  double *r = grid.r;

  a_r_lapl[0] = c_r_lapl[grid.Nr - 1] = 0.;

  for (int i = 1; i < grid.Nr; ++i) {
	a_r_lapl[i] = (3. - r[i + 1] / r[i]) / ((r[i + 1] - r[i - 1]) * (r[i] - r[i - 1]));
  }

  for (int i = 0; i < grid.Nr; ++i) {
	b_r_lapl[i] = -a_r_lapl[i] - c_r_lapl[i];
  }

  for (int i = 0; i < grid.Nr - 1; ++i) {
	c_r_lapl[i] = (3. - r[i - 1] / r[i]) / ((r[i + 1] - r[i - 1]) * (r[i + 1] - r[i]));
  }
}

void set_lower_diag_lapl(size_t num, complex double *low_lapl, double *a_lapl, double coefficient) {
  for (int i = 0; i < num; ++i) {
	low_lapl[i] = coefficient * a_lapl[i];
  }
}

void set_upper_diag_lapl(size_t num, complex double *up_lapl, double *c_lapl, double coefficient) {
  for (int i = 0; i < num; ++i) {
	up_lapl[i] = coefficient * c_lapl[i];
  }
}

double *a_t_lapl;
double *b_t_lapl;
double *c_t_lapl;

double *a_r_lapl;
double *b_r_lapl;
double *c_r_lapl;

void create_abc_arrays(size_t size, double **a, double **b, double **c) {
  *a = create_array(size * sizeof(double));
  *b = create_array(size * sizeof(double));
  *c = create_array(size * sizeof(double));
};

void set_t_laplacian_coefficients(struct grid_t_r_z grid, double *a_t_lapl, double *b_t_lapl, double *c_t_lapl) {
  double dt = grid.t[0] - grid.t[-1];
  for (int it = 0; it < grid.Nt - 1; ++it) {
	a_t_lapl[it] = 1. / square(dt);
  }

  for (int it = 0; it < grid.Nt; ++it) {
	b_t_lapl[it] = -2. / square(dt);
  }

  for (int it = 0; it < grid.Nt - 1; ++it) {
	c_t_lapl[it] = 1. / square(dt);
  }
}

void set_diagonal_lapl(double complex *diagonal_lapl, const double *a_lapl, const double *b_lapl, const double *c_lapl,
					   double complex const_lapl, double dz, int num, double coefficient) {
  for (int i = 0; i <= num - 1; ++i) {
	diagonal_lapl[i] = coefficient * b_lapl[i] + (const_lapl / dz);
  }
  diagonal_lapl[0] += coefficient * a_lapl[0];
  diagonal_lapl[num - 1] += coefficient * c_lapl[num - 1];
}

void tridiagonal_solve(int size, double complex *l, double complex *d, double complex *u, double complex **x,
					   double complex *v) {
  Utmp[0] = u[0] / u[0];

  *x[0] = v[0] / d[0];

  for (int i = 1; i <= size - 1; ++i) {
	double complex tmp = (d[i] - l[i] * Utmp[i - 1]);
	Utmp[i] = u[i] / tmp;
	*x[i] = (v[i] - l[i] * *x[i - 1]) / tmp;
  }

  for (int i = size - 2; i >= 0; --i) {
	*x[i] -= Utmp[i] * *x[i + 1];
  }
}