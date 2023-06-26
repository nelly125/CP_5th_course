#include <iostream>
#include <ctime>
#include "./src/godunov_solver.hpp"

#include "Tests/Piston_functions.hpp"
#include "Tests/course_paper_tests.hpp"

void solve_parameters(double amplitude,
                      double omega, double sigma,
                      double l, double time, uint32_t N,
                      const std::function<double(double x_0,
                                                 double left_0,
                                                 double amplitude,
                                                 double omega,
                                                 double time)> &P_l_func,
                      const std::function<double(double x_0,
                                                 double left_0,
                                                 double amplitude,
                                                 double omega,
                                                 double time)> &P_r_func,
                      boundaries left_boundary,
                      boundaries right_boundary) {

  double x_0 = 2;
  double x_n = 3;
  double x_c = x_0 + l;
  double r_l = 1;
  double u_l = 0;
  double p_l = 1;
  double r_r = r_l * sigma;
  double u_r = 0;
  double p_r = 1;

  std::string dir_name;
  std::ofstream time_log;
  time_log.open("./results/time.log", std::ofstream::out | std::ofstream::app);
  if (!time_log.is_open()) {
    throw std::runtime_error("can't open file");
  }

  time_log << std::scientific;
  double start_time;

  start_time = clock();
  dir_name = solver::solve_system(x_0, x_n, x_c, time, N,
                                  r_l, u_l, p_l, r_r, u_r, p_r,
                                  amplitude, omega, P_l_func, P_r_func,
                                  left_boundary, right_boundary, false);
  time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}

void simple_test(double amplitude,
                 double omega,
                 double l, double time, uint32_t N,
                 const std::function<double(double x_0,
                                            double left_0,
                                            double amplitude,
                                            double omega,
                                            double time)> &P_l_func,
                 const std::function<double(double x_0,
                                            double left_0,
                                            double amplitude,
                                            double omega,
                                            double time)> &P_r_func,
                 boundaries left_boundary,
                 boundaries right_boundary) {

  double x_0 = 0;
  double x_n = 1;

  double x_c = x_0 + l;
  double r_l = 1;
  double u_l = 0;
  double p_l = 1;
  double r_r = r_l;
  double u_r = u_l;
  double p_r = p_l;

/*  double x_c = 0.3;
  double r_l = 1;
  double u_l = 0;
  double p_l = 1;
  double r_r = 1;
  double u_r = 0;
  double p_r = 1;*/

/*  double x_c = 0.5;
  double r_l = 1;
  double u_l = 0;
  double p_l = 1;
  double r_r = 0.125;
  double u_r = 0;
  double p_r = 0.1;*/

  std::string dir_name;
  std::ofstream time_log;
  time_log.open("./results/time.log", std::ofstream::out | std::ofstream::app);
  if (!time_log.
          is_open()
          ) {
    throw std::runtime_error("can't open file");
  }

  time_log <<
           std::scientific;
  double start_time;

  start_time = clock();
/*  dir_name   = solver::solve_system(x_0, x_n, x_c, time, N,
                                    r_l, u_l, p_l, r_r, u_r, p_r,
                                    amplitude, omega, P_l_func, P_r_func,
                                    left_boundary, right_boundary, false);*/

  dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                             r_l, u_l, p_l,
                                             amplitude, omega, P_l_func, P_r_func,
                                             left_boundary, right_boundary, true, limiter::none);

  time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}


double P_1_5_constant(double x_0, double left_0, double amplitude, double omega, double time) {
  return 1.5;
}

double P_1_2_constant(double x_0, double left_0, double amplitude, double omega, double time) {
  return 1.2;
}

double P_0_1_constant(double x_0, double left_0, double amplitude, double omega, double time) {
  return 0.1;
}

double U_0_5_constant(double x_0, double left_0, double amplitude, double omega, double time) {
  return 0.5;
}



void two_piston_wave() {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1;
  double u_l = 0;
  double p_l = 1. / GAMMA;

  double amplitude = 0.01;
  double omega = 2 * M_PI;

  double time = 70;

  auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                  r_l, u_l, p_l,
                                                  amplitude, omega, P_wave_function, P_wave_function,
                                                  boundaries::piston, boundaries::piston, false,
                                                  limiter::none);
}

void piston_up() {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1;
  double u_l = 0;
  double p_l = 1. / GAMMA;

/*
  double amplitude = 0.01;
  double omega = 2 * M_PI;
*/

  double amplitude = 0.01;
  double omega = 50;

  double time = 0.2;

  auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                  r_l, u_l, p_l,
                                                  amplitude, omega, U_0_3_constant, P_1_GAMMA_constant,
                                                  boundaries::piston_U, boundaries::hard_wall, false,
                                                  limiter::none);
}

void piston_down() {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1;
  double u_l = 0;
  double p_l = 1. / GAMMA;

  double amplitude = 0.01;
  double omega = 2 * M_PI;

  double time = 0.5;

  auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                  r_l, u_l, p_l,
                                                  amplitude, omega, P_0_1_constant, P_wave_function,
                                                  boundaries::piston, boundaries::soft, false,
                                                  limiter::none);
}

void two_waves_minmod() {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1;
  double u_l = 0;
  double p_l = 1. / GAMMA;

  double amplitude = 0.2;
  double omega = 2 * M_PI * 5;

  double time = 0.5;

  auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                  r_l, u_l, p_l,
                                                  amplitude, omega, P_function_1_GAMMA_one_positive_wave,
                                                  P_wave_function,
                                                  boundaries::piston, boundaries::piston, false,
                                                  limiter::none);
}


void left_waves() {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1;
  double u_l = 0;
  double p_l = 1. / GAMMA;

  double amplitude = 0.2;
  double omega = 50;

  double time = 100;

  auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                  r_l, u_l, p_l,
                                                  amplitude, omega, P_wave_function, P_1_GAMMA_constant,
                                                  boundaries::piston, boundaries::piston, false,
                                                  limiter::none);
}

void sod_test() {
  double x_0 = 0;
  double x_n = 1;
  uint32_t N = 1000;

  double r_l = 1.6923076923076923;
  double u_l = 0.6067798762169179;
  double p_l = 1.5;

  double r_r = 1;
  double u_r = 0;
  double p_r = 1/GAMMA;

  double x_c = 0.5;


  double amplitude = 0.2;
  double omega = 2 * M_PI * 5;

  double time = 0.2;

  auto dir_name = solver::solve_system(x_0, x_n, x_c, time, N,
                                       r_l, u_l, p_l,
                                       r_r, u_r, p_r,
                                       amplitude, omega, {}, {},
                                       boundaries::soft, boundaries::soft, true);
}

int main() {

/*
  solve_parameters(0.3,
//                   2 * M_PI * sqrt(GAMMA) / 0.3 ,
                  50,
                   25,
                   0.3,
                   150,
                   1000,
                   P_function,
                   P_constant,
                   boundaries::piston,
                   boundaries::piston);
*/

/*    simple_test(0.2,
//              2 * M_PI * 50,
                2 * M_PI * 5,
                0.3,
                0.6,
                1000,
                P_constant,
                P_constant,
                boundaries::piston,
                boundaries::piston);*/

//    piston_wave();
//  piston_down();
//  two_waves_minmod();

//  left_waves();
//  sod_test();


//  piston_up_P();
//  piston_up_U();
//  left_waves_P();
//  left_waves_U(0.03,  M_PI * 2);
    left_waves_P(0.001,  M_PI * 2, 1000);

  return 0;


}
