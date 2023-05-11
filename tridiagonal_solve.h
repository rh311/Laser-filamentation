#ifndef FILAMENTATION_REFACTORING_TRIDIAGONAL_SOLVE_H
#define FILAMENTATION_REFACTORING_TRIDIAGONAL_SOLVE_H

extern double *a_t_lapl;
extern double *b_t_lapl;
extern double *c_t_lapl;

extern double *a_r_lapl;
extern double *b_r_lapl;
extern double *c_r_lapl;

void create_abc_arrays(size_t size, double **a, double **b, double **c);

void tridiagonal_solve_init(int max_matrix_size);

void tridiagonal_solve(int size,
					   double complex *L,
					   double complex *D,
					   double complex *U,
					   double complex **X,
					   double complex *B);

void set_r_laplacian_coefficients(struct grid_t_r_z grid, double *a_r_lapl, double *b_r_lapl, double *c_r_lapl);

void set_t_laplacian_coefficients(struct grid_t_r_z grid, double *a_t_lapl, double *b_t_lapl, double *c_t_lapl);
void set_diagonal_lapl(double complex *diagonal_t_lapl,
					   const double *a_t_lapl,
					   const double *b_t_lapl,
					   const double *c_t_lapl,
					   double complex const_t_lapl,
					   double dz,
					   int Nt,
					   double abc_multiplier);

void set_lower_diag_lapl(size_t num, complex double *low_lapl, double *a_lapl, double coefficient);
void set_upper_diag_lapl(size_t num, complex double *up_t_lapl, double *c_t_lapl, double coefficient);

#endif
