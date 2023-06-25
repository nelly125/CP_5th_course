#ifndef GODUNOV_SOLVER_CONSTANTS_HPP
#define GODUNOV_SOLVER_CONSTANTS_HPP

#define GAMMA (5.0/3.0)
#define EPS		1.0e-10
#define MAX_ITER	20
#define CFL		0.2
#define EPS_CONV 1e-9

const double
        g1 = (GAMMA - 1.0)/(2.0*GAMMA),
        g2 = (GAMMA + 1.0)/(2.0*GAMMA),
        g3 = 2.0*GAMMA/(GAMMA - 1.0),
        g4 = 2.0/(GAMMA - 1.0),
        g5 = 2.0/(GAMMA + 1.0),
        g6 = (GAMMA - 1.0)/(GAMMA + 1.0),
        g7 = (GAMMA - 1.0)/2.0,
        g8 = GAMMA - 1.0;


#endif //GODUNOV_SOLVER_CONSTANTS_HPP
