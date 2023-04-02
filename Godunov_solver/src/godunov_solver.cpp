#include <cmath>
#include <iostream>
#include <cassert>

#include "godunov_solver.hpp"
#include "Riemann_problem/starpu.hpp"
#include "system/system_helper.hpp"
#include "Riemann_problem/sample.hpp"
#include "HLLC/hllc.hpp"

double minmod(double a, double b) {
  if (a * b <= 0.0) {
    return 0.0;
  } else if (std::abs(a) < std::abs(b)) {
    return a;
  } else {
    return b;
  }
}

double minmod_limiter(double phi_im1,
                      double phi_i,
                      double phi_ip1,
                      double dx_im1,
                      double dx_i,
                      double dx_ip1,
                      double u_im1,
                      double u_i,
                      double u_ip1) {
  double r1 = (phi_i - phi_im1) / dx_im1 - (u_i - u_im1) / 2.0;
  double r2 = (phi_ip1 - phi_i) / dx_i - (u_ip1 - u_i) / 2.0;
  double theta = minmod(r1, r2);

  double a = phi_i - phi_im1;
  double b = phi_ip1 - phi_i;
  double sm = 0.0;
  if (a * b > 0.0) {
    double am = std::abs(a);
    double bm = std::abs(b);
    double cm = std::abs(theta * dx_i);
    if (am < bm && am < cm) {
      sm = std::copysign(am, a);
    } else if (bm < am && bm < cm) {
      sm = std::copysign(bm, b);
    } else {
      sm = theta * dx_i;
    }
  }

  return phi_i + 0.5 * sm;
}

double vanleer_limiter(double phi_im1,
                       double phi_i,
                       double phi_ip1,
                       double dx_im1,
                       double dx_i,
                       double dx_ip1,
                       double u_im1,
                       double u_i,
                       double u_ip1) {
  double r1 = (phi_i - phi_im1) / dx_im1 - (u_i - u_im1) / 2.0;
  double r2 = (phi_ip1 - phi_i) / dx_i - (u_ip1 - u_i) / 2.0;
  double theta = 0.0;
  if (r1 * r2 > 0.0) {
    theta = (2.0 * r1 * r2) / (r1 + r2);
  }
  return phi_i + 0.5 * dx_i * theta;
}



double superbee_limiter(double phi_im1,
                        double phi_i,
                        double phi_ip1,
                        double dx_im1,
                        double dx_i,
                        double dx_ip1,
                        double u_im1,
                        double u_i,
                        double u_ip1) {
    double r1 = 2.0 * (phi_i - phi_im1) / (dx_i + dx_im1) - (u_i - u_im1);
    double r2 = 2.0 * (phi_ip1 - phi_i) / (dx_i + dx_ip1) - (u_ip1 - u_i);
    double s1 = (std::abs(r1) < std::abs(r2)) ? r1 : r2;
    double s2 = 2.0 * (phi_ip1 - phi_im1) / (dx_i + dx_ip1 + dx_im1);
    double s3 = 0.0;
    if (s1 * s2 >= 0.0) {
        s3 = s1;
    } else if (std::abs(s1) > std::abs(s2)) {
        s3 = s2;
    } else {
        s3 = s1;
    }
    return phi_i + 0.5 * dx_i * s3;
}

double van_albada(double a, double b) {
  double s = a * b;
  if (s <= 0.0) {
    return 0.0;
  } else {
    return s / (a*a + b*b);
  }
}

double van_albada_limiter(double phi_im1,
                          double phi_i,
                          double phi_ip1,
                          double dx_im1,
                          double dx_i,
                          double dx_ip1,
                          double u_im1,
                          double u_i,
                          double u_ip1) {
  double r1 = (phi_i - phi_im1) / dx_im1 - (u_i - u_im1) / 2.0;
  double r2 = (phi_ip1 - phi_i) / dx_i - (u_ip1 - u_i) / 2.0;
  double theta = van_albada(r1, r2);

  double a = phi_i - phi_im1;
  double b = phi_ip1 - phi_i;
  double sm = 0.0;
  if (a * b > 0.0) {
    double am = std::abs(a);
    double bm = std::abs(b);
    double cm = std::abs(theta * dx_i);
    if (am < bm && am < cm) {
      sm = std::copysign(am, a);
    } else if (bm < am && bm < cm) {
      sm = std::copysign(bm, b);
    } else {
      sm = theta * dx_i;
    }
  }

  return phi_i + 0.5 * sm;
}

void gas_parameters_minmod(uint32_t k,
                           std::vector<gas_parameters> &phi,
                           std::vector<double> &dx,
                           std::vector<double> &u,
                           double &r_int, double &u_int, double &p_int) {
    uint32_t N = phi.size();

    double dx_im1 = k > 0 ? dx[k - 1] : dx[k];
    double dx_i = dx[k];
    double dx_ip1 = k < N - 1 ? dx[k + 1] : dx[k];

    double u_im1 = k > 0 ? u[k - 1] : 0;
    double u_i = u[k];
    double u_ip1 = k < N - 1 ? u[k + 1] : 0;

    {
        double phi_im1 = k > 0 ? phi[k - 1].r : phi[k].r;
        double phi_i = phi[k].r;
        double phi_ip1 = k < N - 1 ? phi[k + 1].r : phi[k].r;
        r_int = van_albada_limiter(phi_im1, phi_i, phi_ip1, dx_im1, dx_i, dx_ip1, u_im1, u_i, u_ip1);
    }

    {
        double phi_im1 = k > 0 ? phi[k - 1].u : phi[k].u;
        double phi_i = phi[k].u;
        double phi_ip1 = k < N - 1 ? phi[k + 1].u : phi[k].u;
        u_int = van_albada_limiter(phi_im1, phi_i, phi_ip1, dx_im1, dx_i, dx_ip1, u_im1, u_i, u_ip1);
    }

    {
        double phi_im1 = k > 0 ? phi[k - 1].p : phi[k].p;
        double phi_i = phi[k].p;
        double phi_ip1 = k < N - 1 ? phi[k + 1].p : phi[k].p;
        p_int = van_albada_limiter(phi_im1, phi_i, phi_ip1, dx_im1, dx_i, dx_ip1, u_im1, u_i, u_ip1);
    }

}

double solver::compute_dt(const std::vector<gas_parameters> &gas,
                          const uint32_t N,
                          const uint32_t i_c,
                          double left_dx,
                          double right_dx) {
    double dx;
    double dt = 1e6;
    for (uint32_t i = 0; i < N - 1; ++i) {
        if (i < i_c) {
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

double solver::compute_dt(const std::vector<gas_parameters> &gas,
                          const uint32_t N,
                          double dx) {
    double dt = 1e6;
    for (uint32_t i = 0; i < N - 1; ++i) {
        double c_l = sqrt(GAMMA * gas[i].p / gas[i].r);
        double c_r = sqrt(GAMMA * gas[i + 1].p / gas[i + 1].r);

        double s_l = fmin(gas[i].u, gas[i + 1].u) - fmax(c_l, c_r);
        double s_r = fmax(gas[i].u, gas[i + 1].u) + fmax(c_l, c_r);
        dt = fmin(dx / fmax(fabs(s_l), fabs(s_r)) * CFL, dt);
    }
    return dt;
}

std::string solver::solve_system(double x_0,
                                 double x_n,
                                 double x_c,
                                 double time,
                                 uint32_t N,
                                 double r_l,
                                 double u_l,
                                 double p_l,
                                 double r_r,
                                 double u_r,
                                 double p_r,
                                 double _amplitude,
                                 double _omega,
                                 const std::function<double(double x_0,
                                                            double left_0,
                                                            double amplitude,
                                                            double omega,
                                                            double time)> &left_piston_func,
                                 const std::function<double(double x_0,
                                                            double left_0,
                                                            double amplitude,
                                                            double omega,
                                                            double time)> &right_piston_func,
                                 boundaries left_boundary,
                                 boundaries right_boundary,
                                 bool hllc_flag) {

    std::vector<gas_parameters> gas_0(N);
    std::vector<gas_parameters> gas_1(N);

    double left_0 = x_0;
    double right_0 = x_n;
    double diaph_0;

    double n_cells = N;
    bool contact = true;

    uint32_t i_contact;
    double s_left;
    double s_right;
    double s_contact;

    double diaph_1;
    double left_1;
    double right_1;

    double dx_left_0, dx_right_0;
    double dx_left_1, dx_right_1;

    double dt_1;

    uint32_t n, m;

    double Energy_0, Energy_1, delta_Energy;
    double A_left_piston, A_right_piston;
    double A_both_pistons;

    double kinetic_energy, internal_energy;

    double amplitude;
    double omega;

    diaph_0 = x_c;
    if (diaph_0 < left_0 || right_0 < left_0) {
        throw std::runtime_error("invalid data");
    }
    i_contact = 0;

    amplitude = _amplitude;
    omega = _omega;

    for (uint32_t i = 0; i < N; ++i) {
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

    n = i_contact;
    m = N - i_contact;

    dx_left_0 = (diaph_0 - left_0) / n;
    dx_right_0 = (right_0 - diaph_0) / m;

    Energy_0 = total_energy(gas_0, i_contact, dx_left_0, dx_right_0);
    Energy_1 = Energy_0;
    delta_Energy = Energy_1 - Energy_0;
    A_left_piston = 0;
    A_right_piston = 0;

    kinetic_energy = total_kinetic_energy(gas_0);
    internal_energy = total_internal_energy(gas_0);

    std::string dir = mk_dir(N, time, amplitude, omega, r_r, x_c);
    std::string data_dir = dir + "/data/";
    std::ofstream output_parameters;
    output_parameters.open(data_dir + "piston_r_u_p.txt");

    std::ofstream cell_params;
    cell_params.open(data_dir + "cell_params.txt");

    std::ofstream output_trajectories;
    output_trajectories.open(data_dir + "trajectories.txt");

    std::ofstream output_delta_energy;
    output_delta_energy.open(data_dir + "delta_energy.txt");

    std::ofstream output_energy;
    output_energy.open(data_dir + "energy.txt");

    double ctime = 0;

    uint32_t step = 0;
    double P_1, P_2, P_3, PP_1, PP_2, PP_3;
    double s_cell;

    double r_left, u_left, p_left;
    double r_right, u_right, p_right;

    bool minmod_flag = false;

    double dx_1, dx_2;
    double xi_0, xi_1;

    double U_0 = 0, P_0, a_0;
    double U_n = 0, P_n, a_n;

    uint32_t output_N = 100;
    if (time >= 10) {
        output_N = N * 10;
    }
    if (N < 1000) {
        output_N = 1;
    }

    double output_time = 0;
    double cell_output_time = 0;

    std::vector<double> s_i(N + 1);
    std::vector<double> x_i_0(N);
    std::vector<double> dx_i_0(N);
    std::vector<double> x_i_1(N);
    std::vector<double> dx_i_1(N);
    std::vector<flows> B(N + 1);

    bool is_contact = false;

    while (ctime < time) {

        if (step % output_N == 0) {
            output_parameters_to_file(output_parameters, gas_0, ctime, left_0, diaph_0, right_0, i_contact, N);

            output_time += 1;

            std::cout << ctime << std::endl;
        }
        if (ctime > cell_output_time) {
            trajectory_to_file(output_trajectories, ctime, left_0, diaph_0, right_0);

            output_cell_params_to_file(cell_params, gas_0, ctime, left_0, diaph_0, right_0, i_contact, N);

            output_delta_energy << ctime << "\t" << delta_Energy << "\t" << A_left_piston << "\t" << A_right_piston
                                << std::endl;
            cell_output_time += 0.05;
        }

        dt_1 = compute_dt(gas_0, N, i_contact, dx_left_0, dx_right_0);
        if (ctime + dt_1 > time) {
            dt_1 = time - ctime;
        }

        if (left_boundary == boundaries::piston) {
            P_0 = left_piston_func(x_0, left_0, amplitude, omega, ctime);
            if (P_0 >= gas_0[0].p) {
                a_0 = sqrt(gas_0[0].r * ((GAMMA + 1) / 2 * P_0 + (GAMMA - 1) / 2 * gas_0[0].p));
                U_0 = gas_0[0].u + (P_0 - gas_0[0].p) / a_0;
            } else {
                U_0 = gas_0[0].u - 2 / (GAMMA - 1) * sqrt(GAMMA * gas_0[0].p / gas_0[0].r) *
                                   (1 - pow(P_0 / gas_0[0].p, g1));
            }
        } else {
            U_0 = 0;
            P_0 = 0;
        }
        if (right_boundary == boundaries::piston) {
            P_n = right_piston_func(x_0, left_0, amplitude, omega, ctime);
            if (P_n >= gas_0[N - 1].p) {
                a_n = sqrt(gas_0[N - 1].r * ((GAMMA + 1) / 2 * P_n + (GAMMA - 1) / 2 * gas_0[N - 1].p));
                U_n = gas_0[N - 1].u - (P_n - gas_0[N - 1].p) / a_n;
            } else {
                U_n = gas_0[N - 1].u + 2 / (GAMMA - 1) * sqrt(GAMMA * gas_0[N - 1].p / gas_0[N - 1].r) *
                                       (1 - pow(P_n / gas_0[N - 1].p, g1));
            }
        } else {
            U_n = 0;
            P_n = 0;
        }

        s_left = U_0;
        s_contact = find_s_cell(gas_0[i_contact], gas_0[i_contact + 1]);
        s_right = U_n;

        left_1 = left_0 + s_left * dt_1;
        diaph_1 = diaph_0 + s_contact * dt_1;
        right_1 = right_0 + s_right * dt_1;

        dx_left_1 = (diaph_1 - left_1) / i_contact;
        dx_right_1 = (right_1 - diaph_1) / (N - i_contact);

        for (uint32_t i = 0; i < N; ++i) {
            x_i_0[i] = find_x(i, left_0, diaph_0, right_0, N, i_contact);
            x_i_1[i] = find_x(i, left_1, diaph_1, right_1, N, i_contact);
            dx_i_0[i] = i < i_contact ? dx_left_0 : dx_right_0;
            dx_i_1[i] = i < i_contact ? dx_left_1 : dx_right_1; //TODO unused
        }

        s_i[0] = U_0;
        s_i[N] = U_n;

        for (uint32_t i = 1; i < N; ++i) {
            s_i[i] = (x_i_1[i] - x_i_0[i]) / dt_1;
        }

        for (uint32_t i = 0; i < N; ++i) {
            if (i == 0) {
                std::cout << i << std::endl;
            }
            if (i == 1) {
                std::cout << i << std::endl;
            }
            if (i == N - 1) {
                std::cout << i << std::endl;
            }

            P_1 = P_2 = P_3 = 0.0;
            PP_1 = PP_2 = PP_3 = 0.0;

            if (i < i_contact) {
                dx_1 = dx_left_0;
                dx_2 = dx_left_1;
            } else {
                dx_1 = dx_right_0;
                dx_2 = dx_right_1;
            }

            if (right_boundary == boundaries::piston && i == N - 1) {
                P_1 = 0;
                P_2 = P_n;
                P_3 = P_n * U_n;
                B[N] = flows(P_1, P_2, P_3);
            } else if (right_boundary == boundaries::soft && i == N - 1) {
                P_1 = 0;
                P_2 = gas_0[N - 1].p;
                P_3 = gas_0[N - 1].p * gas_0[N - 1].u;
            } else if (right_boundary == boundaries::hard_wall && i == N - 1) {
                P_1 = 0;
                P_2 = 0;
                P_3 = 0;
            } else {
                xi_0 = find_x(i + 1, left_0, diaph_0, right_0, N, i_contact);
                xi_1 = find_x(i + 1, left_1, diaph_1, right_1, N, i_contact);

                assert(xi_0 >= left_0 || xi_0 <= right_0);
                assert(xi_1 >= left_1 && xi_1 <= right_1);

                s_cell = s_i[i + 1];

                if (i + 1 == i_contact) {
                    assert(fabs(s_cell - s_contact) < 1e-7);
                }

                r_left = gas_0[i].r;
                u_left = gas_0[i].u;
                p_left = gas_0[i].p;

                r_right = gas_0[i + 1].r;
                u_right = gas_0[i + 1].u;
                p_right = gas_0[i + 1].p;

/*        gas_parameters_minmod(i, gas_0, dx_i_0, s_i, r_left, u_left, p_left);
        gas_parameters_minmod(i + 1, gas_0, dx_i_0, s_i, r_right, u_right, p_right);*/

                if (hllc_flag) {
                    hllc(r_left, u_left, p_left, r_right, u_right, p_right, dx_1, 1.0, s_cell, P_1, P_2, P_3);
                } else {
                    sample(s_cell, r_left, u_left, p_left,
                           r_right, u_right, p_right,
                           1.0, dx_1, P_1, P_2, P_3);
                }
            }

            PP_1 = P_1;
            PP_2 = P_2;
            PP_3 = P_3;

            P_1 = P_2 = P_3 = 0.0;

            if (left_boundary == boundaries::piston && i == 0) {
                P_1 = 0;
                P_2 = P_0;
                P_3 = P_0 * U_0;
            } else if (left_boundary == boundaries::soft && i == 0) {
                P_1 = 0;
                P_2 = gas_0[0].p;
                P_3 = gas_0[0].p * gas_0[0].u;
            } else if (left_boundary == boundaries::hard_wall && i == 0) {
                P_1 = 0;
                P_2 = 0;
                P_3 = 0;
            } else {
                xi_0 = find_x(i, left_0, diaph_0, right_0, N, i_contact);
                xi_1 = find_x(i, left_1, diaph_1, right_1, N, i_contact);

                assert(xi_0 >= left_0 || xi_0 <= right_0);
                assert(xi_1 >= left_1 && xi_1 <= right_1);

                s_cell = s_i[i];

                if (i == i_contact) {
                    assert(fabs(s_cell - s_contact) < 1e-7);
                }

                r_left = gas_0[i - 1].r;
                u_left = gas_0[i - 1].u;
                p_left = gas_0[i - 1].p;

                r_right = gas_0[i].r;
                u_right = gas_0[i].u;
                p_right = gas_0[i].p;
/*        gas_parameters_minmod(i - 1, gas_0, dx_i_0, s_i, r_left, u_left, p_left);
        gas_parameters_minmod(i, gas_0, dx_i_0, s_i, r_right, u_right, p_right);*/

                if (hllc_flag) {
                    hllc(r_left, u_left, p_left, r_right, u_right, p_right, dx_1, 1.0, s_cell, P_1, P_2, P_3);
                } else {
                    sample(s_cell, r_left, u_left, p_left,
                           r_right, u_right, p_right,
                           1.0, dx_1, P_1, P_2, P_3);
                }
            }

            B[i] = flows(P_1, P_2, P_3);

            PP_1 -= P_1;
            PP_2 -= P_2;
            PP_3 -= P_3;

            gas_1[i].r = gas_0[i].r * dx_1 / dx_2 - dt_1 * PP_1 / dx_2;
            gas_1[i].u = (gas_0[i].r * gas_0[i].u * dx_1 / dx_2 - dt_1 * PP_2 / dx_2) / gas_1[i].r;
            gas_1[i].p = ((-dt_1 * PP_3 / dx_2
                           + (gas_0[i].p / (GAMMA - 1.0) + 0.5 * gas_0[i].r * gas_0[i].u * gas_0[i].u) * dx_1 / dx_2) -
                          0.5 * gas_1[i].r * gas_1[i].u * gas_1[i].u) * (GAMMA - 1.0);
        }

        //TODO
        if (right_boundary == boundaries::soft) {
            gas_1[N - 1] = gas_1[N - 2];
        } else if (right_boundary == boundaries::hard_wall) {
            gas_1[N - 1] = gas_1[N - 2];
            gas_1[N - 1].u = -gas_1[N - 1].u;
        }

        for (uint32_t i = 0; i < N; ++i) {
            gas_0[i].r = gas_1[i].r;
            gas_0[i].u = gas_1[i].u;
            gas_0[i].p = gas_1[i].p;
        }

        Energy_1 = total_energy(gas_0, i_contact, dx_left_1, dx_right_1);
        kinetic_energy = total_kinetic_energy(gas_0);
        internal_energy = total_internal_energy(gas_0);

        if (ctime > 0 && left_boundary == boundaries::piston && right_boundary == boundaries::piston) {
            delta_Energy = (Energy_1 - Energy_0) / dt_1;
            A_left_piston = P_0 * U_0;
            A_right_piston = P_n * U_n;
            A_both_pistons = A_left_piston - A_right_piston;
            assert(fabs(delta_Energy - A_both_pistons) < 1e-6);
        }

        diaph_0 = diaph_1;
        left_0 = left_1;
        right_0 = right_1;
        dx_left_0 = dx_left_1;
        dx_right_0 = dx_right_1;
        Energy_0 = Energy_1;

        ctime += dt_1;
        std::cout << ctime << std::endl;

        step++;
    }

    output_parameters_to_file(output_parameters, gas_0, ctime, left_0, diaph_0, right_0, i_contact, N);
    trajectory_to_file(output_trajectories, ctime, left_0, diaph_0, right_0);
    output_delta_energy << ctime << "\t" << delta_Energy << "\t" << A_left_piston << "\t" << A_right_piston
                        << std::endl;
    output_energy << ctime << "\t" << kinetic_energy << "\t" << internal_energy << std::endl;
    std::cout << ctime << std::endl;

    plots(dir, N);

    return dir;
}

double solver::find_s_cell(const gas_parameters left, const gas_parameters right) {
    double p_c;
    double u_c;
    starpu(p_c, u_c, left.r, left.u, left.p, right.r, right.u, right.p);
    return u_c;
}

double solver::find_x(uint32_t j,
                      const double left,
                      const double diaph,
                      const double right,
                      uint32_t N,
                      uint32_t i_contact) {
    if (j < i_contact) {
        return left + double(j) / i_contact * (diaph - left);
    } else {
        return diaph + double(j - i_contact) / (N - i_contact) * (right - diaph);
    }
}

void solver::output_parameters_to_file(std::ostream &out,
                                       const std::vector<gas_parameters> &gas,
                                       double time,
                                       double left,
                                       double diaph,
                                       double right,
                                       uint32_t i_contact,
                                       uint32_t N) {
    out << std::scientific;
    double xpos;

    for (uint32_t i = 0; i < N; ++i) {
        if (i < i_contact) {
            xpos = left + double(i) / i_contact * (diaph - left);
        } else {
            xpos = diaph + double(i - i_contact) / (N - i_contact) * (right - diaph);
        }
        out << xpos << "\t" << gas[i] << time;
        out << std::endl;
    }
}

void solver::output_parameters_to_file(std::ostream &out,
                                       const std::vector<gas_parameters> &gas,
                                       double time,
                                       double left,
                                       double dx,
                                       uint32_t N) {
    out << std::scientific;
    double xpos;

    for (uint32_t i = 0; i < N; ++i) {
        xpos = left + i * dx;
        out << xpos << "\t" << gas[i] << time;
        out << std::endl;
    }
}

void solver::output_cell_params_to_file(std::ostream &out,
                                        const std::vector<gas_parameters> &gas,
                                        double time,
                                        double left,
                                        double diaph,
                                        double right,
                                        uint32_t i_contact,
                                        uint32_t N) {
    out << std::scientific;
    out << time << "\t" << gas[i_contact] << std::endl;
}

void solver::trajectory_to_file(std::ostream &out,
                                double time,
                                double left,
                                double diaph,
                                double right) {
    out << std::scientific;
    out << time << "\t" << left << "\t" << diaph << "\t" << right << std::endl;
}

void solver::trajectory_to_file(std::ostream &out,
                                double time,
                                double left,
                                double right) {
    out << std::scientific;
    out << time << "\t" << left << "\t" << "\t" << right << std::endl;
}

double solver::total_energy(const std::vector<gas_parameters> &gas,
                            uint32_t i_contact,
                            double dx_left,
                            double dx_right) {
    double dx;
    double energy = 0;
    for (uint32_t i = 0; i < gas.size(); ++i) {
        if (i < i_contact) {
            dx = dx_left;
        } else {
            dx = dx_right;
        }
        energy += gas[i].get_total_energy() * dx;
    }
    return energy;
}

double solver::total_energy(const std::vector<gas_parameters> &gas,
                            double dx) {
    double energy = 0;
    for (const auto &ga: gas) {
        energy += ga.get_total_energy() * dx;
    }
    return energy;
}

double solver::total_kinetic_energy(const std::vector<gas_parameters> &gas) {
    double kinetic_energy = 0;
    for (const auto &gas_i: gas) {
        kinetic_energy += gas_i.get_kinetic_energy();
    }
    return kinetic_energy;
}

double solver::total_internal_energy(const std::vector<gas_parameters> &gas) {
    double internal_energy = 0;
    for (const auto &gas_i: gas) {
        internal_energy += gas_i.get_internal_energy();
    }
    return internal_energy;
}

std::string solver::solve_system_boundaries(double x_0,
                                            double x_n,
                                            double time,
                                            uint32_t N,
                                            double r_,
                                            double u_,
                                            double p_,
                                            double _amplitude,
                                            double _omega,
                                            const std::function<double(double x_0,
                                                                       double left_0,
                                                                       double amplitude,
                                                                       double omega,
                                                                       double time)> &left_piston_func,
                                            const std::function<double(double x_0,
                                                                       double left_0,
                                                                       double amplitude,
                                                                       double omega,
                                                                       double time)> &right_piston_func,
                                            boundaries left_boundary,
                                            boundaries right_boundary,
                                            bool hllc_flag) {

    std::vector<gas_parameters> gas_0(N);
    std::vector<gas_parameters> gas_1(N);

    double left_0 = x_0;
    double right_0 = x_n;

    double n_cells = N;
    bool contact = true;

    double s_left;
    double s_right;

    double left_1;
    double right_1;

    double dx_0;
    double dx_1;

    double dt_1;

    double Energy_0, Energy_1, delta_Energy;
    double A_left_piston, A_right_piston;
    double A_both_pistons;

    double kinetic_energy, internal_energy;

    double amplitude;
    double omega;

    if (right_0 < left_0) {
        throw std::runtime_error("invalid data");
    }

    amplitude = _amplitude;
    omega = _omega;

    std::vector<double> s_i(N + 1);
    std::vector<double> x_i_0(N + 1);
    std::vector<double> dx_i_0(N);
    std::vector<double> x_i_1(N + 1);
    std::vector<double> dx_i_1(N);
    std::vector<flows> B(N + 1);
    std::vector<cons_flows> u_0(N);
    std::vector<cons_flows> u_1(N);

    for (uint32_t i = 0; i < N; ++i) {
        double dx = (right_0 - left_0) / n_cells;
        double xpos = left_0 + i * dx + 0.5 * dx;
        gas_0[i].r = r_;
        gas_0[i].u = u_;
        gas_0[i].p = p_;

        u_0[i] = gas_to_cons_flow(gas_0[i]);
    }

    dx_0 = (right_0 - left_0) / N;

    Energy_0 = total_energy(gas_0, dx_0);
    Energy_1 = Energy_0;
    delta_Energy = Energy_1 - Energy_0;
    A_left_piston = 0;
    A_right_piston = 0;

    kinetic_energy = total_kinetic_energy(gas_0);
    internal_energy = total_internal_energy(gas_0);

    std::string dir = mk_dir(N, time, amplitude, omega, r_, 0.);
    std::string data_dir = dir + "/data/";
    std::ofstream output_parameters;
    output_parameters.open(data_dir + "piston_r_u_p.txt");

    std::ofstream cell_params;
    cell_params.open(data_dir + "cell_params.txt");

    std::ofstream output_trajectories;
    output_trajectories.open(data_dir + "trajectories.txt");

    std::ofstream output_delta_energy;
    output_delta_energy.open(data_dir + "delta_energy.txt");

    std::ofstream output_energy;
    output_energy.open(data_dir + "energy.txt");

    double ctime = 0;

    uint32_t step = 0;
    double P_1, P_2, P_3, PP_1, PP_2, PP_3;
    double s_cell;

    double r_left, u_left, p_left;
    double r_right, u_right, p_right;

    bool minmod_flag = false;

    double xi_0, xi_1;

    double U_0 = 0, P_0, a_0;
    double U_n = 0, P_n, a_n;

    uint32_t output_N = 100;
    if (time >= 10) {
        output_N = N * 10;
    }
    if (N < 1000) {
        output_N = 1;
    }

    double output_time = 0;
    double cell_output_time = 0;

    bool is_contact = false;

    while (ctime < time) {

        if (step % output_N == 0) {
            output_parameters_to_file(output_parameters, gas_0, ctime, left_0, dx_0, N);

            output_time += 1;

            std::cout << ctime << std::endl;
        }
        if (ctime > cell_output_time) {
            trajectory_to_file(output_trajectories, ctime, left_0, right_0);

            output_delta_energy << ctime << "\t" << delta_Energy << "\t" << A_left_piston << "\t" << A_right_piston
                                << std::endl;
            cell_output_time += 0.05;
        }

        dt_1 = compute_dt(gas_0, N, dx_0);
        if (ctime + dt_1 > time) {
            dt_1 = time - ctime;
        }

        if (left_boundary == boundaries::piston) {
            P_0 = left_piston_func(x_0, left_0, amplitude, omega, ctime);
            if (P_0 >= gas_0[0].p) {
                a_0 = sqrt(gas_0[0].r * ((GAMMA + 1) / 2 * P_0 + (GAMMA - 1) / 2 * gas_0[0].p));
                U_0 = gas_0[0].u + (P_0 - gas_0[0].p) / a_0;
            } else {
                U_0 = gas_0[0].u - 2 / (GAMMA - 1) * sqrt(GAMMA * gas_0[0].p / gas_0[0].r) *
                                   (1 - pow(P_0 / gas_0[0].p, g1));
            }
        } else {
            U_0 = 0;
            P_0 = 0;
        }
        if (right_boundary == boundaries::piston) {
            P_n = right_piston_func(x_0, left_0, amplitude, omega, ctime);
            if (P_n >= gas_0[N - 1].p) {
                a_n = sqrt(gas_0[N - 1].r * ((GAMMA + 1) / 2 * P_n + (GAMMA - 1) / 2 * gas_0[N - 1].p));
                U_n = gas_0[N - 1].u - (P_n - gas_0[N - 1].p) / a_n;
            } else {
                U_n = gas_0[N - 1].u + 2 / (GAMMA - 1) * sqrt(GAMMA * gas_0[N - 1].p / gas_0[N - 1].r) *
                                       (1 - pow(P_n / gas_0[N - 1].p, g1));
            }
        } else {
            U_n = 0;
            P_n = 0;
        }

        s_left = U_0;
        s_right = U_n;

        left_1 = left_0 + s_left * dt_1;
        right_1 = right_0 + s_right * dt_1;

        dx_1 = (right_1 - left_1) / N;

        for (uint32_t i = 0; i < N + 1; ++i) {
            x_i_0[i] = left_0 + i * dx_0;
            x_i_1[i] = left_1 + i * dx_1;
            if (i == N) {
                continue;
            }
            dx_i_0[i] = dx_0;
            dx_i_1[i] = dx_1; //TODO unused
        }

        s_i[0] = U_0;
        s_i[N] = U_n;
        for (uint32_t i = 1; i < N; ++i) {
            s_i[i] = (x_i_1[i] - x_i_0[i]) / dt_1;
        }

        if (right_boundary == boundaries::piston) {
            P_1 = 0;
            P_2 = P_n;
            P_3 = P_n * U_n;
        } else if (right_boundary == boundaries::soft) {
            P_1 = 0;
            P_2 = gas_0[N - 1].p;
            P_3 = gas_0[N - 1].p * gas_0[N - 1].u;
        } else if (right_boundary == boundaries::hard_wall) {
            P_1 = 0;
            P_2 = 0;
            P_3 = 0;
        }
        B[N] = flows(P_1, P_2, P_3);


        if (left_boundary == boundaries::piston) {
            P_1 = 0;
            P_2 = P_0;
            P_3 = P_0 * U_0;
        } else if (left_boundary == boundaries::soft) {
            P_1 = 0;
            P_2 = gas_0[0].p;
            P_3 = gas_0[0].p * gas_0[0].u;
        } else if (left_boundary == boundaries::hard_wall) {
            P_1 = 0;
            P_2 = 0;
            P_3 = 0;
        }
        B[0] = flows(P_1, P_2, P_3);


        for (uint32_t i = 1; i < N; ++i) {
            if (i == 1) {
              std::cout << "";
            }

          if (i == N - 1) {
            std::cout << "";
          }

            xi_0 = left_0 + i * dx_0;
            xi_1 = left_1 + i * dx_1;

            assert(xi_0 >= left_0 || xi_0 <= right_0);
            assert(xi_1 >= left_1 && xi_1 <= right_1);

            s_cell = s_i[i];

/*            r_left = gas_0[i - 1].r;
            u_left = gas_0[i - 1].u;
            p_left = gas_0[i - 1].p;
            r_right = gas_0[i].r;
            u_right = gas_0[i].u;
            p_right = gas_0[i].p;*/

            gas_parameters_minmod(i - 1, gas_0, dx_i_0, s_i, r_left, u_left, p_left);
            gas_parameters_minmod(i, gas_0, dx_i_0, s_i, r_right, u_right, p_right);

            if (hllc_flag) {
                hllc(r_left, u_left, p_left, r_right, u_right, p_right, dx_1, 1.0, s_cell, P_1, P_2, P_3);
            } else {
                sample(s_cell, r_left, u_left, p_left,
                       r_right, u_right, p_right,
                       1.0, dx_1, P_1, P_2, P_3);
            }

            B[i] = flows(P_1, P_2, P_3);

        }
        for (uint32_t i = 0; i < N; ++i) {
            double temp = dx_0 / dx_1;
            u_1[i].m = (u_0[i].m * temp) - ((dt_1 * (B[i + 1].M - B[i].M))) / (dx_1);
            u_1[i].j = (u_0[i].j * temp) - ((dt_1 * (B[i + 1].J - B[i].J))) / (dx_1);
            u_1[i].e = (u_0[i].e * temp) - ((dt_1 * (B[i + 1].E - B[i].E))) / (dx_1);

            gas_1[i] = c_flow_to_gas(u_1[i]);
        }


        //TODO
        if (right_boundary == boundaries::soft) {
            gas_1[N - 1] = gas_1[N - 2];
        } else if (right_boundary == boundaries::hard_wall) {
            gas_1[N - 1] = gas_1[N - 2];
            gas_1[N - 1].u = -gas_1[N - 1].u;
        }

        for (uint32_t i = 0; i < N; ++i) {
            gas_0[i].r = gas_1[i].r;
            gas_0[i].u = gas_1[i].u;
            gas_0[i].p = gas_1[i].p;
        }

        Energy_1 = total_energy(gas_1, dx_1);
        kinetic_energy = total_kinetic_energy(gas_0);
        internal_energy = total_internal_energy(gas_0);

        if (ctime > 0 && left_boundary == boundaries::piston && right_boundary == boundaries::piston) {
            delta_Energy = (Energy_1 - Energy_0) / dt_1;
            A_left_piston = P_0 * U_0;
            A_right_piston = P_n * U_n;
            A_both_pistons = A_left_piston - A_right_piston;
            assert(fabs(delta_Energy - A_both_pistons) < 1e-6);
        }

        left_0 = left_1;
        right_0 = right_1;
        dx_0 = dx_1;
        Energy_0 = Energy_1;
        u_0 = u_1;

        ctime += dt_1;

        step++;
    }

    output_parameters_to_file(output_parameters, gas_0, ctime, left_0, dx_0, N);
    trajectory_to_file(output_trajectories, ctime, left_0, right_0);
    output_delta_energy << ctime << "\t" << delta_Energy << "\t" << A_left_piston << "\t" << A_right_piston
                        << std::endl;
    output_energy << ctime << "\t" << kinetic_energy << "\t" << internal_energy << std::endl;
    std::cout << ctime << std::endl;

    plots(dir, N);

    return dir;
}
