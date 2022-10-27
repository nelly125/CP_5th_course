#ifndef GODUNOV_SOLVER_GAS_PARAMETERS_HPP
#define GODUNOV_SOLVER_GAS_PARAMETERS_HPP

#include <cmath>
#include "./constants.hpp"

struct gas_parameters {
  double r;
  double u;
  double p;

  double get_total_energy() const {
    return pow(u, 2.) + p / ((GAMMA - 1) * r);
  }

  double get_kinetic_energy() const {
    return pow(u, 2.);
  }

  double get_internal_energy() const {
    return p / ((GAMMA - 1) * r);
  }

  friend std::ostream &operator<<(std::ostream &os, const gas_parameters &gas);
};


#endif //GODUNOV_SOLVER_GAS_PARAMETERS_HPP
