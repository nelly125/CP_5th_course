#include <cmath>
#include <iostream>
#include <cassert>
#include "godunov_solver.hpp"
#include "Riemann_problem/starpu.hpp"
#include "system/system_helper.hpp"
#include "Riemann_problem/sample.hpp"
#define sgn(x) (x>0 ? 1. : -1.)

/*void solver::solve_system(double x_0, double x_n, double x_c, double time, int N, double r, double u, double p) {

}*/

double minmod(double x, double y) {
  if (x * y > 0) {
    return sgn(x) * fmin(fabs(x), fabs(y));
  } else {
    return 0;
  }
}

double solver::compute_dt(const std::vector<gas_parameters> &gas,
                          const int N,
                          const int i_c,
                          double left_dx,
                          double right_dx) {
  double dx;
  double dt = 1e6;
  for (int i = 0; i < N - 1; ++i) {
    if (i <= i_c) {
      dx = left_dx;
    } else {
      dx = right_dx;
    }
    double c_l = sqrt(GAMMA * gas[i].p / gas[i].r);
    double c_r = sqrt(GAMMA * gas[i + 1].p / gas[i + 1].r);

    double s_l = fmin(gas[i].u, gas[i + 1].u) - fmax(c_l, c_r);
    double s_r = fmax(gas[i].u, gas[i + 1].u) + fmax(c_l, c_r);
    dt = fmin(dx / fmax(fabs(s_l), fabs(s_r)) * CFL, dt);
  }
  return dt;
}

void solver::solve_system(double x_0, double x_n, double x_c, double t, int N, double r_l, double u_l, double p_l,
                          double r_r, double u_r, double p_r) {

  std::vector<gas_parameters> gas_0(N);
  std::vector<gas_parameters> gas_1(N);

  double time = t;
  double left_0 = x_0;
  double right_0 = x_n;
  double n_cells = N;
  double diaph_0;
  bool contact = true;

  int i_contact;
  double s_left;
  double s_right;
  double s_contact;

  double diaph_1;
  double left_1;
  double right_1;

  double dx_left_0, dx_right_0;
  double dx_left_1, dx_right_1;

  double dt_1;

  int n, m;

  if (x_c > 0) {
    diaph_0 = x_c;
  }

  for (int i = 1; i < N; ++i) {
    double dx = (right_0 - left_0) / n_cells;
    double xpos = left_0 + i * dx + 0.5 * dx;
    if (xpos <= diaph_0) {
      gas_0[i].r = r_l;
      gas_0[i].u = u_l;
      gas_0[i].p = p_l;
    } else {
      if (contact) {
        i_contact = i;
        contact = false;
      }
      gas_0[i].r = r_r;
      gas_0[i].u = u_r;
      gas_0[i].p = p_r;
    }
  }
  contact = true;

  gas_0[0].r = r_l;
  gas_0[0].u = 0.0;
  gas_0[0].p = p_l;

  gas_0[N - 1].r = r_r;
  gas_0[N - 1].u = 0.0;
  gas_0[N - 1].p = p_r;

  n = i_contact;
  m = N - i_contact;

  dx_left_0 = (diaph_0 - left_0) / n;
  dx_right_0 = (right_0 - diaph_0) / m;

  std::string dir = mk_dir(N, time);
  std::ofstream output_parameters;
  output_parameters.open(dir + "piston_r_u_p.txt");

  double ctime = 0;

  int step = 0;
  double P_1, P_2, P_3, PP_1, PP_2, PP_3;
  double s_cell;

  double r_left, u_left, p_left;
  double r_right, u_right, p_right;

  bool minmod_flag = false;

  double dx_1, dx_2;
  double xi_0, xi_1;

  while (ctime < time) {
    if (step % 100 == 0) {
      output_parameters_to_file(output_parameters, gas_0, ctime, left_0, diaph_0, right_0, i_contact, N);
      std::cout << ctime << std::endl;
    }

    dt_1 = compute_dt(gas_0, N, i_contact, dx_left_0, dx_right_0);

    gas_0[0].p = 1 + 0.5 * sin(100* (ctime + dt_1));
//    gas_0[0].p = 0.5;


    s_left = find_s_cell(gas_0[0], gas_0[1]);
    s_contact = find_s_cell(gas_0[i_contact - 1], gas_0[i_contact]);
    s_right = find_s_cell(gas_0[N - 2], gas_0[N - 1]);

    left_1 = left_0 + s_left * dt_1;
    diaph_1 = diaph_0 + s_contact * dt_1;
    right_1 = right_0 + s_right * dt_1;


    dx_left_1 = (diaph_1 - left_1) / n;
    dx_right_1 = (right_1 - diaph_1) / m;

    for (int i = 1; i < N - 1; ++i) {
      P_1 = P_2 = P_3 = 0.0;
      PP_1 = PP_2 = PP_3 = 0.0;

      if (i <= i_contact) {
        dx_1 = dx_left_0;
        dx_2 = dx_left_1;
      } else {
        dx_1 = dx_right_0;
        dx_2 = dx_right_1;
      }

      xi_0 = find_x(i + 1, left_0, diaph_0, right_0, N, i_contact);
      xi_1 = find_x(i + 1, left_1, diaph_1, right_1, N, i_contact);

      assert(xi_0 >= left_0 || xi_0 <= right_0);
      assert(xi_1 >= left_1 && xi_1 <= right_1);

      s_cell = (xi_1 - xi_0) / dt_1;

      if (i == N) {
        assert(fabs(s_cell - s_right) < 1e-7);
      }

      get_values(i, N, minmod_flag, gas_0[i], r_left, u_left, p_left);
      get_values(i + 1, N, minmod_flag, gas_0[i + 1], r_right, u_right, p_right);
      sample(s_cell, r_left, u_left, p_left,
             r_right, u_right, p_right,
             1.0, dx_1, P_1, P_2, P_3);

      PP_1 = P_1;
      PP_2 = P_2;
      PP_3 = P_3;

      P_1 = P_2 = P_3 = 0.0;

      xi_0 = find_x(i, left_0, diaph_0, right_0, N, i_contact);
      xi_1 = find_x(i, left_1, diaph_1, right_1, N, i_contact);
      assert(xi_0 >= left_0 || xi_0 <= right_0);
      assert(xi_1 >= left_1 && xi_1 <= right_1);
      s_cell = (xi_1 - xi_0) / dt_1;
      if (i == 1) {
        assert(fabs(s_cell - s_left) < 1e-7);
      }

      get_values(i - 1, N, minmod_flag, gas_0[i - 1], r_left, u_left, p_left);
      get_values(i, N, minmod_flag, gas_0[i], r_right, u_right, p_right);

      sample(s_cell, r_left, u_left, p_left,
             r_right, u_right, p_right,
             1.0, dx_1, P_1, P_2, P_3);

      PP_1 -= P_1;
      PP_2 -= P_2;
      PP_3 -= P_3;

      gas_1[i].r = gas_0[i].r * dx_1 / dx_2 - dt_1 * PP_1 / dx_2;
      gas_1[i].u = (gas_0[i].u * gas_0[i].r * dx_1 / dx_2 - dt_1 * PP_2 / dx_2) / gas_1[i].r;
      gas_1[i].p = ((-dt_1 * PP_3 / dx_2 +
          (gas_0[i].p / (GAMMA - 1.0) + 0.5 * gas_0[i].r * gas_0[i].u * gas_0[i].u) * dx_1 / dx_2) -
          0.5 * gas_1[i].r * gas_1[i].u * gas_1[i].u) * (GAMMA - 1.0);
    }

    for (int i = 1; i < N - 1; ++i) {
      gas_0[i].r = gas_1[i].r;
      gas_0[i].u = gas_1[i].u;
      gas_0[i].p = gas_1[i].p;
    }

//    gas_0[0] = gas_0[1];
    gas_0[N - 1] = gas_0[N - 2];

    diaph_0 = diaph_1;
    left_0 = left_1;
    right_0 = right_1;
    dx_left_0 = dx_left_1;
    dx_right_0 = dx_right_1;

    if (ctime + dt_1 > time) {
      dt_1 = time - ctime;
    }
    ctime += dt_1;
//    std::cout << ctime << std::endl;
    step++;

  }

  plots(dir, N);
}

double solver::find_s_cell(const gas_parameters left, const gas_parameters right) {
  double p_c;
  double u_c;
  starpu(p_c, u_c, left.r, left.u, left.p, right.r, right.u, right.p);
  return u_c;
}

double solver::find_x(int j, const double left, const double diaph, const double right, int N, int i_contact) {
  if (j < 1) {
    throw std::runtime_error("");
  }
  if (j <= i_contact) {
    return left + double(j - 1) / i_contact * (diaph - left);
  } else {
    return diaph + double(j - 1 - i_contact) / (N - i_contact) * (right - diaph);
  }
}

void solver::get_values(int k,
                        int N,
                        bool minmod_flag,
                        gas_parameters gas,
                        double &r,
                        double &u,
                        double &p) {
  if (k == 0 || k == N - 1) {
    r = gas.r;
    u = gas.u;
    p = gas.p;
  } else {
    r = gas.r + int(minmod_flag) * 1. / 2 * minmod(gas.r - gas.r, gas.r - gas.r);
    u = gas.u + int(minmod_flag) * 1. / 2 * minmod(gas.u - gas.u, gas.u - gas.u);
    p = gas.p + int(minmod_flag) * 1. / 2 * minmod(gas.p - gas.p, gas.p - gas.p);
  }

}
void solver::output_parameters_to_file(std::ostream &out,
                                       const std::vector<gas_parameters> gas,
                                       double time,
                                       double left,
                                       double diaph,
                                       double right,
                                       int i_contact,
                                       int N) {
  double xpos;
  for (int i = 0; i < N; ++i) {
    if (i <= i_contact) {
      xpos = left + double(i) / i_contact * (diaph - left);
    } else {
      xpos = diaph + double(i - i_contact) / (N - i_contact) * (right - diaph);
    }
    out << xpos << "\t" << gas[i];
    if (i == 0) {
      out << time;
    }
    out << std::endl;
  }
}

