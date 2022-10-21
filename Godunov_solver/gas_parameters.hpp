#ifndef GODUNOV_SOLVER_GAS_PARAMETERS_HPP
#define GODUNOV_SOLVER_GAS_PARAMETERS_HPP

#include <cmath>

struct gas_parameters {
    double r;
    double u;
    double p;

    double get_total_energy() {
        return pow(u, 2.) + p / ((GAMMA - 1) * r);
    }
};

#endif //GODUNOV_SOLVER_GAS_PARAMETERS_HPP
