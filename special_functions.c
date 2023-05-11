#include <math.h>

#define NMAX 6
#define H 0.4
#define A1 (2.0/3.0)
#define A2 0.4
#define A3 (2.0/7.0)

double h = 1.05457e-27;
extern double pi;
extern double electronMass;
extern double electronCharge;
extern double omega;

double dawson(double x) {
  int i, n0;
  double d1, d2, e1, e2, sum, x2, xp, xx, ans, tmp;
  static double c[NMAX + 1];
  static int init = 0;
  if (init == 0) {
	init = 1;
	for (i = 1; i <= NMAX; i++) {
	  tmp = (2. * i - 1.) * H;
	  tmp *= tmp;
	  c[i] = exp(-tmp);
	}
  }
  if (fabs(x) < 0.2) {
	x2 = x * x;
	ans = x * (1.0 - A1 * x2 * (1.0 - A2 * x2 * (1.0 - A3 * x2)));
  } else {
	xx = fabs(x);
	n0 = 2 * (int)(0.5 * xx / H + 0.5);
	xp = xx - n0 * H;
	e1 = exp(2.0 * xp * H);
	e2 = e1 * e1;
	d1 = n0 + 1;
	d2 = d1 - 2.0;
	sum = 0.0;
	for (i = 1; i <= NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
	  sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));

	ans = 0.5641895835 * exp(-xp * xp) * sum;
	if (x < 0.)
	  ans = -ans;
  }
  return ans;
}

double EllipticK(double y) {
  double ERRTOL = 0.08;
  double C1 = (1.0 / 24.0);
  double C2 = 0.1;
  double C3 = (3.0 / 44.0);
  double C4 = (1.0 / 14.0);

  double alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt, zt;

  xt = 0.;
  yt = 1. - y;
  zt = 1.;
  do {
	sqrtx = sqrt(xt);
	sqrty = sqrt(yt);
	sqrtz = sqrt(zt);
	alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
	xt = 0.25 * (xt + alamb);
	yt = 0.25 * (yt + alamb);
	zt = 0.25 * (zt + alamb);
	ave = (xt + yt + zt) / 3.0;
	delx = (ave - xt) / ave;
	dely = (ave - yt) / ave;
	delz = (ave - zt) / ave;
  } while (fmax(fmax(fabs(delx), fabs(dely)), fabs(delz)) > ERRTOL);
  e2 = delx * dely - delz * delz;
  e3 = delx * dely * delz;
  return (1.0 + (C1 * e2 - C2 - C3 * e3) * e2 + C4 * e3) / sqrt(ave);
}

double rd(double y) {
  double ERRTOL = 0.05;
  double C1 = (3.0 / 14.0);
  double C2 = (1.0 / 6.0);
  double C3 = (9.0 / 22.0);
  double C4 = (3.0 / 26.0);
  double C5 = (0.25 * C3);
  double C6 = (1.5 * C4);

  double alamb, ave, delx, dely, delz, ea, eb, ec, ed, ee, fac, sqrtx, sqrty, sqrtz, sum, xt, yt, zt;

  xt = 0.;
  yt = 1. - y;
  zt = 1.;
  sum = 0.0;
  fac = y / 3.;

  do {
	sqrtx = sqrt(xt);
	sqrty = sqrt(yt);
	sqrtz = sqrt(zt);
	alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
	sum += fac / (sqrtz * (zt + alamb));
	fac = 0.25 * fac;
	xt = 0.25 * (xt + alamb);
	yt = 0.25 * (yt + alamb);
	zt = 0.25 * (zt + alamb);
	ave = 0.2 * (xt + yt + 3.0 * zt);
	delx = (ave - xt) / ave;
	dely = (ave - yt) / ave;
	delz = (ave - zt) / ave;
  } while (fmax(fmax(fabs(delx), fabs(dely)), fabs(delz)) > ERRTOL);

  ea = delx * dely;
  eb = delz * delz;
  ec = ea - eb;
  ed = ea - 6.0 * eb;
  ee = ed + ec + ec;

  return 3. * sum
	  + fac * (1.0 + ed * (-C1 + C5 * ed - C6 * delz * ee) + delz * (C2 * ee + delz * (-C3 * ec + delz * C4 * ea)))
		  / (ave * sqrt(ave));
}

double EllipticE(double k) {
  return EllipticK(k) - rd(k);
}

double Ww(double CnormA, double band_gap) {
  if (CnormA < 7.58926e+9) {
	return 0.;
  }

  double gamm2 = pow(omega / electronCharge, 2) * electronMass * band_gap / CnormA;
  double Eps = 1. / (1. + gamm2);
  double Gamm = (gamm2) / (1. + gamm2);

  double EllipticK_Gamm = EllipticK(Gamm);
  double EllipticE_Gamm = EllipticK_Gamm - rd(Gamm);

  double EllipticK_Eps = EllipticK(Eps);
  double EllipticE_Eps = EllipticK_Eps - rd(Eps);

  double alpha = pi * (EllipticK_Gamm - EllipticE_Gamm) / EllipticE_Eps;
  double beta = pow(pi, 2) / (2. * EllipticK_Eps * EllipticE_Eps);
  double x = 2. * band_gap * EllipticE_Eps / (pi * h * omega * sqrt(Gamm));
  double nu1 = x + 0.5;
  double nu2 = nu1 - x;

  double sum = 0.;
  double intermediate_sum;
  int j = 4;
  int i = 0;

  do {
	intermediate_sum = sum;
	while (i < j) {
	  sum += exp(-i * alpha) * dawson(sqrt(beta * ((double)i + nu2)));
	  ++i;
	}
	j *= 2;
  } while ((intermediate_sum / sum) < 0.99);

  double Q = sqrt(pi / (2. * EllipticK_Eps)) * sum;

  return 2. * (2. * omega / (9. * pi)) * (pow(omega * electronMass / (h * sqrt(Gamm)), 1.5)) * Q * exp(-alpha * nu1);
}