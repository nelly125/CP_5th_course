#include "sample.hpp"
#include "starpu.hpp"
#include "../../constants.hpp"

double sample( double &s, double r_l, double u_l, double p_l, double r_r, double u_r, double p_r, double n, double dx,
               double &F_1, double &F_2, double &F_3, bool flux_flag, bool piston_flag ) {

  double c, cml, cmr, p_cl, p_cr, shl, shr, sl, sr, stl, str;
  double c_l, c_r;
  double d, u, p;
  double p_c, u_c;
  double s_l, s_r;
  double dt;
  double un_l, un_r;

  c_l = sqrt(GAMMA * p_l / r_l);
  c_r = sqrt(GAMMA * p_r / r_r);

  un_l = u_l * n;
  un_r = u_r * n;

  starpu(p_c, u_c, r_l, un_l, p_l, r_r, un_r, p_r);

  if (s <= u_c) {
    if (p_c <= p_l) {
      shl = un_l - c_l;

      if (s <= shl) {
        d = r_l;
        u = un_l;
        p = p_l;
      } else {
        cml = c_l * pow(p_c / p_l, g1);
        stl = u_c - cml;

        if (s > stl) {
          d = r_l * pow(p_c / p_l, 1.0 / GAMMA);
          u = u_c;
          p = p_c;
        } else {
          u = g5 * (c_l + g7 * un_l + s);
          c = g5 * (c_l + g7 * (un_l - s));
          d = r_l * pow(c / c_l, g4);
          p = p_l * pow(c / c_l, g3);
        }
      }
    } else {
      p_cl = p_c / p_l;
      sl = un_l - c_l * sqrt(g2 * p_cl + g1);

      if (s <= sl) {
        d = r_l;
        u = un_l;
        p = p_l;
      } else {
        d = r_l * (p_cl + g6) / (p_cl * g6 + 1.0);
        u = u_c;
        p = p_c;
      }
    }
  } else {
    if (p_c > p_r) {
      p_cr = p_c / p_r;
      sr = un_r + c_r * sqrt(g2 * p_cr + g1);

      if (s >= sr) {
        d = r_r;
        u = un_r;
        p = p_r;
      } else {
        d = r_r * (p_cr + g6) / (p_cr * g6 + 1.0);
        u = u_c;
        p = p_c;
      }
    } else {
      shr = un_r + c_r;

      if (s >= shr) {
        d = r_r;
        u = un_r;
        p = p_r;
      } else {
        cmr = c_r * pow(p_c / p_r, g1);
        str = u_c + cmr;

        if (s <= str) {
          d = r_r * pow(p_c / p_r, 1.0 / GAMMA);
          u = u_c;
          p = p_c;
        } else {
          u = g5 * (-c_r + g7 * un_r + s);
          c = g5 * (c_r - g7 * (un_r - s));
          d = r_r * pow(c / c_r, g4);
          p = p_r * pow(c / c_r, g3);
        }
      }
    }
  }

  if (flux_flag) {
    F_1 = d * (u - s);
    F_2 = F_1 * u + p;
    F_3 = F_1 * (p / ((GAMMA - 1) * d) + u * u / 2) + p * u;
  } else {
    F_1 = d;
    F_2 = u;
    F_3 = p;
  }

  if (piston_flag) {
    u = s;
    F_1 = d * (u - s);
    F_2 = F_1 * u + p;
    F_3 = F_1 * (p / ((GAMMA - 1) * d) + u * u / 2) + p * u;
  }

  s_l = fmin(un_l, un_r) - fmax(c_l, c_r);
  s_r = fmax(un_l, un_r) + fmax(c_l, c_r);
  dt = dx / fmax(fabs(s_l), fabs(s_r)) * CFL;

  return dt;

}


