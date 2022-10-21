#ifndef GODUNOV_SAMPLE_H
#define GODUNOV_SAMPLE_H

#include <cmath>

double sample(double &s, double r_l, double u_l, double p_l, double r_r, double u_r, double p_r, double n, double dx,
              double &F_1, double &F_2, double &F_3, bool flux_flag = true, bool piston_flag = false);


#endif //GODUNOV_SAMPLE_H
