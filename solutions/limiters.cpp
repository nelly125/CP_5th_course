#include <cmath>
#include "limiters.h"
#include "../src/godunov_solver.hpp"

#include <gtest/gtest.h>

namespace limiters {
  double P_1_GAMMA_constant( double x_0, double left_0, double amplitude, double omega, double time ) {
    return 1. / GAMMA;
  }

  double P_1_constant( double x_0, double left_0, double amplitude, double omega, double time ) {
    return 1;
  }

  double P_function_one_wave( double x_0, double left_0, double amplitude, double omega, double time ) {
    double P_0;
    double T = 2 * M_PI / omega;
    if (time < T) {
      P_0 = (1 / GAMMA + amplitude * sin(omega * (time)));
    } else {
      P_0 = 1 / GAMMA;
    }
    return P_0;
  }

  double P_wave_function( double x_0, double left_0, double amplitude, double omega, double time ) {
    return (1 / GAMMA + amplitude * sin(omega * (time)));
  }

  TEST(lLIMITERS, minmod) {

    double   x_0 = 0;
    double   x_n = 1;
    uint32_t N   = 1000;

    double r_l = 1;
    double u_l = 0;
    double p_l = 1. / GAMMA;

    double amplitude = 0.2;
    double omega     = 2 * M_PI * 5;

    double time = 0.7;

    auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                    r_l, u_l, p_l,
                                                    amplitude, omega, P_function_one_wave, P_function_one_wave,
                                                    boundaries::piston, boundaries::piston, false,
                                                    limiter::van_albada_lim);
  }

  TEST(lLIMITERS, one_wave_2_pi) {

    double   x_0 = 0;
    double   x_n = 1;
    uint32_t N   = 1000;

    double r_l = 1;
    double u_l = 0;
    double p_l = 1. / GAMMA;

    double amplitude = 0.2;
    double omega     = 2 * M_PI;

    double time = 10;

    auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                    r_l, u_l, p_l,
                                                    amplitude, omega, P_wave_function, P_1_GAMMA_constant,
                                                    boundaries::piston, boundaries::piston, false,
                                                    limiter::none);
  }

}

int main() {
  RUN_ALL_TESTS();

  return 0;
}