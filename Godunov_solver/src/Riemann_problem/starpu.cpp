#include "starpu.hpp"
#include "../../constants.hpp"

void prefun( double &f, double &fd, double &p,
             double &dk, double &pk, double &ck ) {
  double ak, bk, pratio, qrt;

  if (p <= pk) {
    // волна разрежения
    pratio = p / pk;
    f = g4 * ck * (pow(pratio, g1) - 1.0);
    fd = (1.0 / (dk * ck)) * pow(pratio, -g2);

  } else {
    // ударная волна
    ak = g5 / dk;
    bk = g6 * pk;
    qrt = sqrt(ak / (bk + p));
    f = (p - pk) * qrt;
    fd = (1.0 - 0.5 * (p - pk) / (bk + p)) * qrt;
  }

}

double guessp( double r_l, double u_l, double p_l, double c_l, double r_r, double u_r, double p_r, double c_r ) {
  double cup, gel, ger,
      pmax, pmin, ppv, pq, pm,
      ptl, ptr,
      qmax, quser, um;

  quser = 2.0;
  cup = 0.25 * (r_l + r_r) * (c_l + c_r);
  ppv = 0.5 * (p_l + p_r) + 0.5 * (u_l - u_r) * cup;
  ppv = ppv > 0.0 ? ppv : 0.0;
  pmin = p_l > p_r ? p_r : p_l;
  pmax = p_l > p_r ? p_l : p_r;
  qmax = pmax / pmin;

  if (qmax <= quser && pmin <= ppv && ppv <= pmax) {
    pm = ppv;
  } else {
    if (ppv < pmin) {
      // две волны разрежения
      pq = pow(p_l / p_r, g1);
      um = (pq * u_l / c_l + u_r / c_r + g4 * (pq - 1.0)) / (pq / c_l + 1.0 / c_r);
      ptl = 1.0 + g7 * (u_l - um) / c_l;
      ptr = 1.0 + g7 * (um - u_r) / c_r;
      pm = 0.5 * (p_l * pow(ptl, g3) + p_r * pow(ptr, g3));
    } else {
      // две ударных волны
      gel = sqrt((g5 / r_l) / (g6 * p_l + ppv));
      ger = sqrt((g5 / r_l) / (g6 * p_r + ppv));
      pm = (gel * p_l + ger * p_r - (u_r - u_l)) / (gel + ger);
    }

  }
  return pm;
}

void starpu( double &p, double &u, double r_l, double u_l, double p_l, double r_r, double u_r, double p_r ) {
  int i;

  double pstart,
      pold,
      udiff,
      change;

  double fl, fld,
      fr, frd;

  double c_l, c_r;

  c_l = sqrt(GAMMA * p_l / r_l);
  c_r = sqrt(GAMMA * p_r / r_r);

  pstart = guessp(r_l, u_l, p_l, c_l, r_r, u_r, p_r, c_r);

  pold = pstart;

  udiff = u_r - u_l;

  change = 10.0 * EPS;

  for (i = 0; i < MAX_ITER && change > EPS; i++) {
    prefun(fl, fld, pold, r_l, p_l, c_l);

    prefun(fr, frd, pold, r_r, p_r, c_r);

    p = pold - (fl + fr + udiff) / (fld + frd);

    change = 2.0 * fabs((p - pold) / (p + pold));

    if (p < 0.0) p = 0.0;

    pold = p;
  }

  u = 0.5 * (u_l + u_r + fr - fl);
}
