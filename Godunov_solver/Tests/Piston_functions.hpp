#ifndef GODUNOV_SOLVER_GODUNOV_SOLVER_TESTS_PISTON_FUNCTIONS_HPP_
#define GODUNOV_SOLVER_GODUNOV_SOLVER_TESTS_PISTON_FUNCTIONS_HPP_

#include "../constants.hpp"
#include <math.h>

double P_function_norm( double x_0, double left_0, double amplitude, double omega, double time ) {
  return pow(x_0 / left_0, 2. * GAMMA) * (1.0 + amplitude * sin(omega * (time)));
//  return(1.0 + amplitude * sin(omega * (time)));
}


double P_function( double x_0, double left_0, double amplitude, double omega, double time ) {
  return (1 + amplitude * sin(omega * (time)));
}

double P_constant( double x_0, double left_0, double amplitude, double omega, double time ) {
  return 5;
}

double P_function_one_wave( double x_0, double left_0, double amplitude, double omega, double time ) {
  double P_0;
  double T = 2 * M_PI / omega;
  if (time < T) {
    P_0 = (1 + amplitude * sin(omega * (time)));
  } else {
    P_0 = 1;
  }
  return P_0;
}

double P_function_one_positive_wave( double x_0, double left_0, double amplitude, double omega, double time ) {
  double P_0;
  double T = 2 * M_PI / omega;
  if (time < T / 2) {
    P_0 = (1 + amplitude * sin(omega * (time)));
  } else {
    P_0 = 1;
  }
  return P_0;
}

#endif //GODUNOV_SOLVER_GODUNOV_SOLVER_TESTS_PISTON_FUNCTIONS_HPP_
