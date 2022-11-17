#include <iostream>
#include <ctime>
#include "./src/godunov_solver.hpp"

void step_dependency() {
  double amplitude;
  double omega;
  double sigma;
  double l;

  amplitude = 0.1;
  omega = 50;
  sigma = 25;
  l = 0.5;

  double x_0 = 2;
  double x_n = 3;
  double x_c = x_0 + l;
  double time = 100;
  uint32_t N;
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

  for (N = 1000; N <= 1000; N *= 10) {
    start_time = clock();
    dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
    time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

  }
}

void amplitude_dependency() {
  double amplitude;
  double omega;
  double sigma;
  double l;

  amplitude = 0.2;
  omega = 30;
  sigma = 5;
  l = 0.5;

  double x_0 = 2;
  double x_n = 3;
  double x_c = x_0 + l;
  double time = 50;
  uint32_t N;
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

  for (N = 1000; N <= 1000; N *= 10) {
    start_time = clock();
    dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
    time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

  }

}

int main() {
//    solver Solver;
/*

  double amplitude;
  double omega;
  double sigma;
  double l;

  amplitude = 0.2;
  omega = 30;
  sigma = 5;
  l = 0.5;

  double x_0 = 2;
  double x_n = 3;
  double x_c = x_0 + l;
  double time = 50;
  uint32_t N = 1000;
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
*/

/*  omega = 10;
  while (omega >= 10 && omega <= 100) {
    start_time = clock();
    dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
    time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
    omega += 10;
  }*/

/*
  start_time = clock();
  dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
  time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
*/

  step_dependency();

  return 0;

}
