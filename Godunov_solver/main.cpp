#include <iostream>

#include "./src/godunov_solver.hpp"

int main() {
    std::cout << "Hello, World!" << std::endl;

    solver Solver;
    double x_0 = 0;
    double x_n = 1;
    double x_c = 0.5;
    double time = 2;
    int N = 1000;
    double r_l = 1;
    double u_l = 0;
    double p_l = 1;
    double r_r = 1;
    double u_r = 0;
    double p_r = 1;

    Solver.solve_system(x_0, x_n, x_c, time, N, r_l, u_l, p_l, r_r, u_r, p_r);

    return 0;


}
