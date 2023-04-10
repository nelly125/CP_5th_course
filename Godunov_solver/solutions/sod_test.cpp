#include <gtest/gtest.h>

#include "../src/godunov_solver.hpp"
#include <cmath>

/*
TEST(lLIMITERS, SodTest) {

    double x_0 = 0;
    double x_n = 1;
    uint32_t N = 1000;

    double r_l = 1;
    double u_l = 0;
    double p_l = 1;

    double r_r = 0.125;
    double u_r = 0;
    double p_r = 0.1;


    double amplitude = 0.2;
    double omega = 2 * M_PI * 5;

    double time = 0.3;

    auto dir_name = solver::solve_system_boundaries(x_0, x_n, time, N,
                                                    r_r, u_r, p_r,
                                                    amplitude, omega, {}, {},
                                                    boundaries::soft, boundaries::soft, true,
                                                    limiter::none);
}*/
