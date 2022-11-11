#ifndef GODUNOV_STARPU_H
#define GODUNOV_STARPU_H

#include <cmath>
double guessp( double r_l, double u_l, double p_l, double c_l, double r_r, double u_r, double p_r, double c_r );
void starpu( double &p, double &u, double r_l, double u_l, double p_l, double r_r, double u_r, double p_r );
void prefun( double &f, double &fd, double &p,
             double &dk, double &pk, double &ck );

#endif //GODUNOV_STARPU_H
