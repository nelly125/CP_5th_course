#include <iostream>
#include <ctime>
#include <vector>
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

  for (N = 100; N <= 1000; N *= 10) {
    start_time = clock();
    dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
    time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

  }
}

void pressure_dependency() {
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
  double time = 100;
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

  start_time = clock();
  dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
  time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}

void amplitude_dependency() {
  double amplitude;
  double omega;
  double sigma;
  double l;

  omega = 50;
  sigma = 25;
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

  amplitude = 0.25;
  start_time = clock();
  dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
  time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

/*  amplitude = 0.01;
  while (amplitude <= 0.09) {
    start_time = clock();
    dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
    time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
    amplitude += 0.01;
  }*/
}

void omega_dependency() {
  double amplitude;
  double omega;
  double sigma;
  double l;

  amplitude = 0.2;
  sigma = 25;
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

  omega = 10;
  while (omega <= 150) {
    start_time = clock();
    dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
    time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
    omega += 10;
  }
}

void sigma_dependency() {
  double amplitude;
  double omega;
  double sigma;
  double l;

  amplitude = 0.2;
  omega = 50;
  l = 0.5;

  double x_0 = 2;
  double x_n = 3;
  double x_c = x_0 + l;
  double time = 50;
  uint32_t N = 1000;
  double r_l = 1;
  double u_l = 0;
  double p_l = 1;
  double r_r;
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

  sigma = 1;
  while (sigma <= 35) {
    r_r = r_l * sigma;
    start_time = clock();
    dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
    time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
    if (fabs(sigma - 1) < 1e-5) {
      sigma = 5;
    } else {
      sigma += 5;
    }
  }
}

void l_dependency() {
  double amplitude;
  double omega;
  double sigma;
  double l;

  amplitude = 0.2;
  omega = 50;
  sigma = 25;

  double x_0 = 2;
  double x_n = 3;
  double x_c;
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

  l = 0.1;
  while (l <= 1) {
    x_c = x_0 + l;
    start_time = clock();
    dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
    time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
    l += 0.1;
  }
}

void grid_search() {

  double x_0 = 2;
  double x_n = 3;
  double x_c;
  double time = 50;
  uint32_t N = 1000;
  double r_l = 1;
  double u_l = 0;
  double p_l = 1;
  double r_r;
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

  std::vector<double> amplitude_grid = {0.05, 0.1, 0.2, 0.25};
  std::vector<double> omega_grid = {10, 30, 50, 70};
  std::vector<double> l_grid = {0.1, 0.3, 0.5};
  std::vector<double> sigma_grid = {1., 5., 15., 25.};

  for (auto amplitude : amplitude_grid) {
    for (auto omega : omega_grid) {
      for (auto l : l_grid) {
        for (auto sigma : sigma_grid) {
          x_c = x_0 + l;
          r_r = r_l * sigma;
          start_time = clock();
          dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
          time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
        }
      }
    }
  }
}

void solve_parameters( double amplitude,
                       double omega, double sigma,
                       double l, double time, uint32_t N ) {

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
  dir_name = solver::solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r, amplitude, omega);
  time_log << dir_name << ": \t" << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;

}

int main() {

//  pressure_dependency();
//  amplitude_dependency();
//  omega_dependency();
//  sigma_dependency();
//  l_dependency();
//  grid_search();
//  solve_parameters(0.05, 10, 15, 0.1, 200, 1000);
  solve_parameters(0.3, 30, 5, 0.5, 5, 1000);
  return 0;

}
