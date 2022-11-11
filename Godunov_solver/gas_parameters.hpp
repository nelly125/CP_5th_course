#ifndef GODUNOV_SOLVER_GAS_PARAMETERS_HPP
#define GODUNOV_SOLVER_GAS_PARAMETERS_HPP

#include <cmath>
#include "./constants.hpp"

struct gas_parameters {
  double r;
  double u;
  double p;

  double get_total_energy() const {
    return r * (p / (r * (GAMMA - 1)) + pow(u, 2) / 2);
  }

  double get_kinetic_energy() const {
    return r * pow(u, 2.) / 2.;
  }

  double get_internal_energy() const {
    return p / ((GAMMA - 1) );
  }

  double impedans() const {
    return r * sqrt(GAMMA * p / r);
  }

  friend std::ostream &operator<<(std::ostream &os, const gas_parameters &gas);
};


#endif //GODUNOV_SOLVER_GAS_PARAMETERS_HPP
