#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "math_functions.h"
#include "physical_expressions.h"
#include "constants.h"
#include "tridiagonal_solve.h"
#include "create_array.h"
#include <string.h>

double omega;

struct grid_t_r_z grid = {600, 800, 3001};

double dt = 0.8e-15;
double dr = 0.05e-4;
void r_initialization(struct grid_t_r_z grid, double dr) {
  grid.r[0] = dr / 2.0;
  for (int ir = 1; ir < 200; ++ir)
	grid.r[ir] = grid.r[ir - 1] + dr;
  for (int ir = 200; ir < grid.Nr; ++ir) {
	dr *= 1.01;
	grid.r[ir] = grid.r[ir - 1] + dr;
  }
}

void t_initialization(struct grid_t_r_z grid, double dt) {
  for (int i = -1; i <= grid.Nt; ++i)
	grid.t[i] = i * dt;
}

int Npulses = 260;

int ind_rt(int ir, int it) {
  return ir + it * grid.Nr;
}

int ind_rz(int ir, int iz) {
  return ir + iz * grid.Nr;
}

double Ww(double CnormA, double band_gap);

void slice(double complex *input_array, double complex **output_array, int start, int step, int size);

void save(double *array, int number_of_elements, char *file_name) {
  FILE *pf = fopen(file_name, "wb");
  fwrite(array, sizeof(*array), number_of_elements, pf);
  fclose(pf);
}

void read(double *array, int number_of_elements, char *file_name) {
  FILE *pf = fopen(file_name, "rb");
  fread(array, sizeof(*array), number_of_elements, pf);
  fclose(pf);
}

void A_initialization(double complex *A,
					  struct grid_t_r_z grid,
					  double rad_r,
					  double rad_t,
					  double k,
					  double focal_length,
					  double pulseEnergy,
					  double const_intensity) {
  int Nr = grid.Nr;
  int Nt = grid.Nt;
  double *r = grid.r;
  double *t = grid.t;

  double amplitude = sqrt
	  (
		  (pulseEnergy / const_intensity) / (pow(pi, 3. / 2.) * rad_t * square(rad_r))
	  );
  for (int it = 0; it < Nt; ++it) {
	double tmp = gaussian(t[it], t[Nt / 2], rad_t);
	for (int ir = 0; ir < Nr; ++ir) {
	  A[ind_rt(ir, it)] = amplitude * tmp * gaussian(r[ir], 0., rad_r) *
		  cexp(I * k * square(r[ir]) / (2. * focal_length));
	}
  }
}

void B_lapl_initialization(double complex **A,
						   double complex *const_vector_lapl,
						   double complex const_lapl,
						   double dz,
						   int num) {
  for (int i = 0; i < num; ++i)
	const_vector_lapl[i] = (const_lapl / dz) * *A[i];
}

int main() {
  double wavelength = 3.2e-4;
  double pulseEnergy = 250.;
  double pulseWidth = 0.0120;
  double pulseDuration = 100.0e-15;
  //вещество
  double band_gap_ex_eV = 12.8;
  double band_gap_eh_eV = 12.8;//14.3;
  double Nmax = 6.1e+22;
  double tau_col = 10.0e-15;
  double lin_refractive_index = n0(wavelength);
  double GVD = GVD_func(wavelength);
  double n2_I = 8.1e-24;
  double focal_length = 2.5e+20;
  double n_modif_max = 1.;
  double dn;

  double Dt = 1.e-6;
  double thermal_conductivity = 1.42e+6;
  double density = 2.635;
  double heat_capacity = 1.562e+7;
  double aa = thermal_conductivity / (heat_capacity * density);

  double band_gap_ex = band_gap_ex_eV * eV;
  double band_gap_eh = band_gap_eh_eV * eV;
  omega = 2.0 * pi / wavelength * speedOfLight;
  double k = omega * lin_refractive_index / speedOfLight;
  double rad_r = 0.5 * pulseWidth / sqrt(log(2.0));
  double rad_t = 0.5 * pulseDuration / sqrt(log(2.0));
  double n2_A = speedOfLight * lin_refractive_index * n2_I / (8 * pi);

  double Energy[] = {0.9935214819990955, 0.9737377879358617, 0.9485317328189717, 1.0062890843369048, 1.0537184711090848,
					 0.979090813798742, 0.9778721587129245, 1.046880712024718, 0.9608271804239494, 1.0554250430970606,
					 1.0241675150582603, 0.9692908561605319, 0.9998027537231308, 0.9752247753133724, 0.9910191356367499,
					 0.9081518216557138, 0.9293074064478657, 0.9614676091773812, 1.126641588467572, 0.9125310406056163,
					 1.0321686549315265, 1.0221260229864129, 0.9994530246204671, 0.9937775598519916, 0.9585082744106274,
					 0.9788184332473168, 1.0496843855216067, 1.00688753510747, 0.9765133593056735, 1.0527212515447977,
					 0.9677229169683761, 0.9597314582620899, 1.0147017129947398, 1.103963127807276, 1.0219723827507423,
					 1.0811860634463553, 1.086284637479757, 1.0461411645397913, 1.0500034400302112, 0.8787630919969263,
					 0.9802740804615587, 1.0440972566470168, 0.9916546490515311, 1.0775782576356914, 0.9367711713316857,
					 1.026455649197035, 1.0244530479779161, 1.008959403917989, 1.012716895678051, 1.0060731992233731,
					 0.949063408900711, 1.015016755255438, 0.9848709337466809, 1.0552973784947581, 1.0110056496043032,
					 0.8991696482962264, 0.9799002535195555, 0.9153219813038396, 1.0541871943768075, 0.977390797464748,
					 1.0173483323780634, 0.9945263639250044, 0.9935815233759625, 0.9959258128792864, 0.9820551380372493,
					 1.0512193300054617, 1.004575796828621, 0.9971379506909123, 1.0306801836765456, 1.0181723789850148,
					 1.0805190614887885, 0.9884545439410435, 1.0100437381596674, 1.082578674422375, 1.0137397908735823,
					 0.9765322367186189, 1.0873110563144612, 1.0389773233513813, 0.9821628739902745, 1.0316334851369078,
					 1.0708796656900121, 1.0537320157872496, 1.0616614447386272, 0.9839680742644761, 0.9046536717752987,
					 0.9729631870700584, 1.0047016415579524, 1.1538765309210772, 0.8859455361300789, 0.9470453172384583,
					 0.9994158133515926, 0.9952227781722724, 0.9471020633438582, 1.010327801766198, 1.0409067878944034,
					 1.0171319640118524, 0.9653796823881725, 1.025772109214805, 1.0036557779636375, 0.9598994659328209,
					 1.006374645452321, 0.9598167891987115, 0.9717556073057247, 0.9485294507689993, 1.056193131626498,
					 0.9849368734754597, 1.031601184225722, 1.0277074074984573, 1.0207475720564023, 0.9921988509030782,
					 0.9467135261144887, 1.0622979257403113, 1.0327733106437071, 1.0129677232270753, 0.9927500030735615,
					 1.0414758311180972, 0.940425844773944, 0.9602984221679977, 0.991231080577464, 1.0190985646946127,
					 1.0026454302079115, 1.0860953301978566, 0.995811282132177, 1.0572717709049355, 1.0119872159344139,
					 0.9809768931033547, 0.998749505698794, 1.0816741786943094, 1.0258572201547183, 1.0293969101581804,
					 1.0840433671385827, 1.071110291814564, 1.0470481533233738, 1.0131579136084516, 0.9709182997974337,
					 1.017856917707731, 0.9274216281436886, 0.954739726411533, 0.885145412922808, 1.016821051016639,
					 0.9071006885270314, 1.0107952587327762, 0.9779589726352057, 1.002903193040468, 0.9998950582217786,
					 0.9910440101133613, 1.048356034608667, 0.9767020698412721, 0.9921777644514614, 1.0482995124129022,
					 1.0720504736833525, 1.0261775545729133, 0.9289435247206799, 1.0226542777196765, 0.9658578068694658,
					 1.0134886562923633, 0.9408915348136205, 1.0628118219561222, 1.0141724717453622, 0.9964711009445419,
					 0.9994857558990496, 1.0023963627679402, 0.9255830071076014, 1.0210698670563907, 1.0085825776955488,
					 0.980852714504141, 0.9973335648871656, 0.9912564197059415, 1.0219425555353727, 0.9420311018478692,
					 0.9413356084729749, 0.9579967133744249, 0.9826011429539425, 1.072107213035648, 1.0980318988231972,
					 0.9525980349265981, 0.9778618000325932, 0.9974602022527056, 1.003480151614327, 0.9953456460133157,
					 0.9390892369093153, 1.0494951748391395, 0.9454597941787994, 1.0366131489961634, 0.9982843080025282,
					 1.0127225427565292, 1.0180391752698599, 0.9737325181118371, 0.9177138136800832, 1.0010938844871018,
					 1.1085731328025425, 1.06441559315813, 0.9426280097904466, 1.0425953092071363, 0.955292893308834,
					 0.9700834282928555, 1.0601298214108466, 0.9959186545550966, 1.0227878759351685,
					 0.9778721587129245};

  // коэффициенты для использования в расчетах
  double const_intensity = speedOfLight * lin_refractive_index / (8 * pi);
  double complex const_r_lapl = -I * 2. * omega * lin_refractive_index / speedOfLight;
  double complex const_t_lapl = -I * 2.;
  double complex const_self_foc = -I * omega * n2_A / speedOfLight;
  double complex const_modif_foc = -I * omega / speedOfLight;
  double complex const_plasm_defoc =
	  I * 2 * pi * electronCharge * electronCharge / (speedOfLight * omega * lin_refractive_index * electronMass);
  double complex
	  const_avalanche = electronCharge * electronCharge / (2 * band_gap_eh * omega * omega * electronMass * tau_col);
  double complex const_nonlin_absorption_ex = -4 * pi * band_gap_ex / (speedOfLight * lin_refractive_index);
  double complex const_nonlin_absorption_eh = -4 * pi * band_gap_eh / (speedOfLight * lin_refractive_index);
  double complex const_avalanche_absorption = -2 * pi * electronCharge * electronCharge
	  / (speedOfLight * lin_refractive_index * omega * omega * electronMass * tau_col);

  char switch_avalanche = 0;
  char switch_plasm_defoc = 1;
  char switch_n_modif = 1;

  double complex *A = create_array(grid.Nt * grid.Nr * sizeof(double complex));

  double complex **A_r_const_view = create_array(grid.Nt * sizeof(double complex *));

  double complex **A_t_const_view = create_array(grid.Nr * sizeof(double complex *));

  double *A2_rt = create_array(grid.Nt * grid.Nr * sizeof(double));

  grid.r = create_array(grid.Nr * sizeof(double));// r[-1] ... r[Nr]
  r_initialization(grid, dr);

  save(grid.r, grid.Nr, "r_coord.bin");

  grid.t = create_array((grid.Nt + 2) * sizeof(double));// t[-1] ... t[Nt]
  grid.t = grid.t + 1; //t[-1]
  t_initialization(grid, dt);

  create_abc_arrays(grid.Nr, &a_r_lapl, &b_r_lapl, &c_r_lapl);
  set_r_laplacian_coefficients(grid, a_r_lapl, b_r_lapl, c_r_lapl);

  double complex *low_r_lapl = create_array(grid.Nr * sizeof(double complex));
  double complex *diagonal_r_lapl = create_array(grid.Nr * sizeof(double complex));
  double complex *up_r_lapl = create_array(grid.Nr * sizeof(double complex));
  double complex *const_vector_r_lapl = create_array(grid.Nr * sizeof(double complex));

  set_lower_diag_lapl(grid.Nr, low_r_lapl, a_r_lapl, 1.);
  set_upper_diag_lapl(grid.Nr, up_r_lapl, c_r_lapl, 1.);

  create_abc_arrays(grid.Nt, &a_t_lapl, &b_t_lapl, &c_t_lapl);
  set_t_laplacian_coefficients(grid, a_t_lapl, b_t_lapl, c_t_lapl);

  double complex *low_t_lapl = create_array(grid.Nt * sizeof(double complex));
  double complex *diagonal_t_lapl = create_array(grid.Nt * sizeof(double complex));
  double complex *up_t_lapl = create_array(grid.Nt * sizeof(double complex));
  double complex *const_vector_t_lapl = create_array(grid.Nt * sizeof(double complex));

  set_lower_diag_lapl(grid.Nt, low_t_lapl, a_t_lapl, -GVD);
  set_upper_diag_lapl(grid.Nt, up_t_lapl, c_t_lapl, -GVD);

  tridiagonal_solve_init(grid.Nr > grid.Nt ? grid.Nr : grid.Nt);
  double *plasma_rz = create_array(grid.Nr * grid.Nz * sizeof(double));
  double *flux_rz = create_array(grid.Nr * grid.Nz * sizeof(double));
  double *max_intensity_rz = create_array(grid.Nr * grid.Nz * sizeof(double));
  double *z_coord = create_array(grid.Nz * sizeof(double));
  double *n_modif_rz = create_array(grid.Nr * grid.Nz * sizeof(double));
  for (int j = 0; j < grid.Nr * grid.Nz; ++j)
	n_modif_rz[j] = 0.0;

  double *dz = create_array((grid.Nz + 1) * sizeof(double));

  int num_pulses_to_save[] = {1, 2, 4, 8, 16, 32, 64, 128, 256};
  int num_pulses_to_save_cur_ind = 0;

  A_initialization(A, grid, rad_r, rad_t, k, focal_length, pulseEnergy, const_intensity);

  double Anorm_max = 0.0;
  for (int i = 0; i < grid.Nr * grid.Nt; ++i) {
	if (cnorm(A[i]) > Anorm_max)
	  Anorm_max = cnorm(A[i]);
  }

  dz[0] = 2.0e-4;

  double Pcr = Power_cr(wavelength, lin_refractive_index, n2_I);

  double P = 0.;
  for (int j = 0; j < grid.Nr - 1; ++j) { P += cnorm(A[j + (grid.Nt / 2) * grid.Nr]) * grid.r[j]
		* (grid.r[j + 1] - grid.r[j]);
  }
  P *= const_intensity * 2.0 * pi;

  printf("time_length/rad_t = %g\n", grid.Nt * dt / rad_t);
  printf("r_max/rad_r = %g\n", grid.r[grid.Nr - 1] / rad_r);
  printf("Pcr = %g\n", Pcr);
  printf("P = %g\n", P);
  printf("P/Pcr = %g\n", P / Pcr);
  printf("z_sf = %g\n", z_sf(omega, lin_refractive_index, rad_r, P / Pcr, focal_length));
  printf("ldiff = %g\n", k * square(rad_r));
  printf("ldisp = %g\n", square(rad_t) / GVD);
  printf("pulse energy = %g\n", pulseEnergy);
  printf("GVD = %g\n", GVD);
  printf("n0 = %g\n", lin_refractive_index);
  printf("const_intensity = %g\n", const_intensity);

  int pulse_count = 0;
  for (int i_pulse = 0; i_pulse < Npulses; ++i_pulse) {
	double z = 0.0;

	A_initialization(A,
					 grid,
					 rad_r,
					 rad_t,
					 k,
					 focal_length,
					 250. * (Energy[i_pulse] * (1. / 5.) + (4. / 5.)) /* pulseEnergy*/ ,
					 const_intensity);

	for (int iz = 0; iz < grid.Nz; ++iz) {
	  //запись A2 в файл
	  if ((iz % 50 == 0)) {
		double max_intens = 0.;
		double CnormA;

		for (int i = 0; i < grid.Nr * grid.Nt; ++i) {
		  CnormA = cnorm(A[i]);
		  A2_rt[i] = CnormA;

		  if (max_intens < CnormA)
			max_intens = CnormA;
		}

		char filename[1024];
		sprintf(filename, "A2_rt_%i_%i.bin", i_pulse, iz);

		save(A2_rt, grid.Nr * grid.Nt, filename);
		printf("step %i\n", iz);
		printf("%lf\n", max_intens);
	  }

	  set_diagonal_lapl(diagonal_r_lapl, a_r_lapl, b_r_lapl, c_r_lapl, const_r_lapl, dz[iz], grid.Nr, 1.);
	  //#pragma omp parallel for schedule(static)
	  for (int it = 0; it < grid.Nt; ++it) {
		slice(A, A_t_const_view, it * grid.Nr, 1, grid.Nr);
		B_lapl_initialization(A_t_const_view, const_vector_r_lapl, const_r_lapl, dz[iz], grid.Nr);
		tridiagonal_solve(grid.Nr, low_r_lapl, diagonal_r_lapl, up_r_lapl, A_t_const_view, const_vector_r_lapl);
	  }

	  for (int ir = 0; ir < grid.Nr; ++ir) {
		double Wtmp;
		double N_eh = 0.;
		double N_ex = 0.;
		//double N_eh_max = 0.;
		double complex w;
		double CnormA;
		for (int it = 0; it < grid.Nt; ++it) {
		  CnormA = cnorm(A[ind_rt(ir, it)]);
		  Wtmp = Ww(CnormA, band_gap_ex);

		  if (Wtmp || N_eh) {
			double tmp = Wtmp;
			if (switch_avalanche)
			  tmp += CnormA * const_avalanche * N_eh;

			double tmp2 = dt * tmp * (1. - (N_eh + N_ex) / Nmax);
			N_eh += 1. * tmp2 * n_modif_rz[ind_rz(ir, iz)];
			N_ex += tmp2 * (1. - 1. * n_modif_rz[ind_rz(ir, iz)]);
		  }
		  w = const_self_foc * CnormA;

		  if (switch_avalanche) {
			w += const_avalanche_absorption * (1.0 - (N_eh /*+ N_ex*/) / Nmax) * N_eh;
		  }

		  if (switch_plasm_defoc) {

			w += const_plasm_defoc * N_eh;
		  }

		  if (switch_n_modif) {
			w += 0.005 * const_modif_foc * n_modif_rz[ir + iz * grid.Nr];
		  }

		  if (CnormA && Wtmp)
			w += const_nonlin_absorption_eh * Wtmp * (1. - (N_eh + N_ex) / Nmax) / CnormA;

		  w *= dz[iz];

		  A[ind_rt(ir, it)] *= cexp(w);
		}

		plasma_rz[ind_rz(ir, iz)] = N_eh + N_ex;

	  }

	  z += dz[iz];

	  if (i_pulse == 0) {
		z_coord[iz] = z;
		dz[iz + 1] = dz[iz];
	  }

	}

	if (i_pulse == 0) {
	  save(z_coord, grid.Nz, "z_coord.bin");
	}

	if (switch_n_modif) {
	  if (i_pulse == 0) {
		double max_Plasma_rz = 0.0;                  //поиск максимального значения концентрации плазмы
		for (int i_zr = 0; i_zr < grid.Nz * grid.Nr; i_zr++)
		  if (plasma_rz[i_zr] > max_Plasma_rz)
			max_Plasma_rz = plasma_rz[i_zr];
		dn = 0.05 * n_modif_max / max_Plasma_rz;
	  }

	  {
		{
		  for (int iz = 0; iz < grid.Nz; ++iz)
			for (int ir = 0; ir < grid.Nr; ++ir)
			  n_modif_rz[ind_rz(ir, iz)] +=
				  plasma_rz[ind_rz(ir, iz)] * dn * (1. - n_modif_rz[ind_rz(ir, iz)] / n_modif_max);
		  ++pulse_count;
		}
	  }
	}

	if (pulse_count == num_pulses_to_save[num_pulses_to_save_cur_ind]) {
	  {
		char filename[1024];
		sprintf(filename, "plasma_rz_%i.bin", pulse_count);
		save(plasma_rz, grid.Nr * grid.Nz, filename);
	  }

	  {
		char filename[1024];
		sprintf(filename, "flux_rz_%i.bin", pulse_count);
		save(flux_rz, grid.Nr * grid.Nz, filename);
	  }
	  //save_array_to_file_one_index(n_modif_rz, grid.Nr * grid.Nz, "n_modif_rz_%i.bin", pulse_count)
	  {
		char filename[1024];
		sprintf(filename, "n_modif_rz_%i.bin", pulse_count);
		save(n_modif_rz, grid.Nr * grid.Nz, filename);
	  }

	  ++num_pulses_to_save_cur_ind;
	}

	double plazma_integral = 0.0;
	for (int iz = 0; iz < grid.Nz; ++iz) {
	  for (int ir = 0; ir < grid.Nr - 1; ++ir) {
		plazma_integral += plasma_rz[ind_rz(ir, iz)] * grid.r[ir] * (grid.r[ir + 1] - grid.r[ir]) * 2.0 * pi * dz[iz];
	  }
	}
	printf("plazma_integral = %g %i \n", plazma_integral, i_pulse);

  }
  free(A);
  free(const_vector_r_lapl);
  free(grid.r);
  free(a_r_lapl);
  free(b_r_lapl);
  free(c_r_lapl);
  free(diagonal_r_lapl);
  free(plasma_rz);
  free(flux_rz);
  free(max_intensity_rz);
  free(dz);
  free(n_modif_rz);
  free(z_coord);
  free(A2_rt);

  return 0;
}

