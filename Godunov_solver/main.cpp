#include <iostream>
#include <ctime>
#include <vector>
#include "./src/godunov_solver.hpp"

#include "Tests/Piston_functions.hpp"

void solve_parameters( double amplitude,
                       double omega, double sigma,
                       double l, double time, uint32_t N,
                       const std::function<double( double x_0,
                                                   double left_0,
                                                   double amplitude,
                                                   double omega,
                                                   double time )> &P_l_func,
                       const std::function<double( double x_0,
                                                   double left_0,
                                                   double amplitude,
                                                   double omega,
                                                   double time )> &P_r_func,
                       boundaries left_boundary,
                       boundaries right_boundary ) {

  double x_0 = 2;
  double x_n = 3;
  double x_c = x_0 + l;
  double r_l = 1;
  double u_l = 0;
  double p_l = 1;
  double r_r = r_l * sigma;
  double u_r = 0;
  double p_r = 1;

  std::string   dir_name;
  std::ofstream time_log;
  time_log.open("./results/time.log", std::ofstream::out | std::ofstream::app);
  if (!time_log.is_open()) {
    throw std::runtime_error("can't open file");
  }

  time_log << std::scientific;
  double start_time;

  start_time = clock();
  dir_name   = solver::solve_system(x_0, x_n, x_c, time, N,
                                    r_l, u_l, p_l, r_r, u_r, p_r,
                                    amplitude, omega, P_l_func, P_r_func,
                                    left_boundary, right_boundary, false);
  time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}

void simple_test( double amplitude,
                  double omega,
                  double l, double time, uint32_t N,
                  const std::function<double( double x_0,
                                              double left_0,
                                              double amplitude,
                                              double omega,
                                              double time )> &P_l_func,
                  const std::function<double( double x_0,
                                              double left_0,
                                              double amplitude,
                                              double omega,
                                              double time )> &P_r_func,
                  boundaries left_boundary,
                  boundaries right_boundary ) {

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

  std::string   dir_name;
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

  dir_name   = solver::solve_system_boundaries(x_0, x_n, time, N,
                                    r_l, u_l, p_l,
                                    amplitude, omega, P_l_func, P_r_func,
                                    left_boundary, right_boundary, true);

  time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

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

  simple_test(0.2,
//              2 * M_PI * 50,
              2 * M_PI*5,
              0.3,
              0.6,
              1000,
              P_constant,
              P_constant,
              boundaries::piston,
              boundaries::piston);

  return 0;

}
