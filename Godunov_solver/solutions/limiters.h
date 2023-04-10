#ifndef GODUNOV_SOLVER_LIMITERS_H
#define GODUNOV_SOLVER_LIMITERS_H
namespace limiters {

    double P_1_GAMMA_constant(double x_0, double left_0, double amplitude, double omega, double time);

    double P_1_constant(double x_0, double left_0, double amplitude, double omega, double time);
}

double P_1_5_constant(double x_0, double left_0, double amplitude, double omega, double time);

double P_function_one_wave(double x_0, double left_0, double amplitude, double omega, double time);

#endif //GODUNOV_SOLVER_LIMITERS_H
