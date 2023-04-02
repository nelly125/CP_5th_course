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
    return p / ((GAMMA - 1));
  }

  double impedans() const {
    return r * sqrt(GAMMA * p / r);
  }

  friend std::ostream &operator<<( std::ostream &os, const gas_parameters &gas );
};

/*
struct cell_flux {
  double m;
  double p;
  double e;
  cell_flux( double _m, double _p, double _e ) : m(_m), p(_p), e(_e) {};
  cell_flux() = default;
};
*/

struct flows {
  flows() = default;
  flows(double d1, double d2, double d3): M(d1), J(d2), E(d3) {};

  double M;
  double J;
  double E;

  flows operator - (flows);
};

struct cons_flows {
  double m;
  double j;
  double e;

  cons_flows( double _m, double _p, double _e ) : m(_m), j(_p), e(_e) {};
  cons_flows() = default;
};


flows gas_to_flows(gas_parameters &gas);

cons_flows gas_to_cons_flow(gas_parameters &gas);

gas_parameters c_flow_to_gas(cons_flows &c_flow);

#endif //GODUNOV_SOLVER_GAS_PARAMETERS_HPP
