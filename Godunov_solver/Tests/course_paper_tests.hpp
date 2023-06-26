#ifndef GODUNOV_SOLVER_COURSE_PAPER_TESTS_HPP
#define GODUNOV_SOLVER_COURSE_PAPER_TESTS_HPP

#include "../constants.hpp"
#include "../src/godunov_solver.hpp"

double P_1_GAMMA_constant(double x_0, double left_0, double amplitude, double omega, double time) {
  return 1. / GAMMA;
}

double U_0_3_constant(double x_0, double left_0, double amplitude, double omega, double time) {
  return 0.32302914123489934;
}

double P_1_constant(double x_0, double left_0, double amplitude, double omega, double time) {
  return 1;
}

void piston_up_U() {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1;
  double u_l = 0;
  double p_l = 1. / GAMMA;

  double amplitude = 0.01;
  double omega = 50;

  double time = 0.2;

  auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                  r_l, u_l, p_l,
                                                  amplitude, omega, U_0_3_constant, P_1_GAMMA_constant,
                                                  boundaries::piston_U, boundaries::soft, false,
                                                  limiter::none);
}

void piston_up_P() {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1;
  double u_l = 0;
  double p_l = 1. / GAMMA;

  double amplitude = 0.01;
  double omega = 50;

  double time = 0.2;

  auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                  r_l, u_l, p_l,
                                                  amplitude, omega, P_1_constant, P_1_GAMMA_constant,
                                                  boundaries::piston_P, boundaries::soft, false,
                                                  limiter::none);
}

double P_wave_function(double x_0, double left_0, double amplitude, double omega, double time) {
  return (1 / GAMMA + amplitude * sin(omega * (time)));
}

double U_wave_function(double x_0, double left_0, double amplitude, double omega, double time) {
  return (amplitude * sin(omega * (time)));
}

void left_waves_P(double amplitude = 0.01, double omega = 2 * M_PI, double time = 200) {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1;
  double u_l = 0;
  double p_l = 1. / GAMMA;

  auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                  r_l, u_l, p_l,
                                                  amplitude, omega, P_wave_function, P_1_GAMMA_constant,
                                                  boundaries::piston_P, boundaries::piston_P, true,
                                                  limiter::none);
}

void left_waves_U(double amplitude = 0.01, double omega = 2 * M_PI, double time = 200) {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1;
  double u_l = 0;
  double p_l = 1. / GAMMA;

  auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                  r_l, u_l, p_l,
                                                  amplitude, omega, U_wave_function, P_1_GAMMA_constant,
                                                  boundaries::piston_U, boundaries::piston_P, true,
                                                  limiter::none);
}



void piston_wave(    double amplitude = 0.01, double omega = 4 * M_PI) {
    double x_0 = 0;
    double x_n = 1;
    uint32_t N = 1000;

    double r_l = 1;
    double u_l = 0;
    double p_l = 1. / GAMMA;



    double time = 200;

    auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                    r_l, u_l, p_l,
                                                    amplitude, omega, P_function_one_wave, P_1_GAMMA_constant,
                                                    boundaries::piston_P, boundaries::piston_P, false,
                                                    limiter::none);
}

#endif //GODUNOV_SOLVER_COURSE_PAPER_TESTS_HPP
