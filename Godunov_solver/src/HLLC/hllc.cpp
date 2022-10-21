#include <cmath>
#include "hllc.hpp"
#include "../../constants.hpp"

double hllc(const double &r_l, const double &u_l, const double &p_l, const double &r_r, const double &u_r,
            const double &p_r, const double &dx, const double &n, const double &s, double &F1, double &F2, double &F3) {
    double s_l, s_r, c_l, c_r;
    double j_l, j_r, D_c, P_c, rho_l, rho_r;
    double un_l, un_r, sn;
    double dt;

    un_l = u_l * n;
    un_r = u_r * n;
    sn = s;

    c_l = sqrt(GAMMA * p_l / r_l);
    c_r = sqrt(GAMMA * p_r / r_r);

    s_l = fmin(un_l, un_r) - fmax(c_l, c_r);
    s_r = fmax(un_l, un_r) + fmax(c_l, c_r);

    j_l = r_l * (un_l - s_l);
    j_r = r_r * (un_r - s_r);
    D_c = (j_r * un_r - j_l * un_l + p_r - p_l) / (j_r - j_l);
    P_c = (j_r * p_l - j_l * p_r - j_l * j_r * (un_r - un_l)) / (j_r - j_l);
    rho_l = r_l * (s_l - un_l) / (s_l - D_c);
    rho_r = r_r * (s_r - un_r) / (s_r - D_c);

    if (s_l >= sn) {
        F1 = r_l * (un_l - sn);
        F2 = (p_l + un_l * F1) * n;
        F3 = F1 * (p_l / ((GAMMA - 1) * r_l) + un_l * un_l / 2) + p_l * un_l;
    } else if (s_l < sn && D_c >= sn) {
        F1 = rho_l * (D_c - sn);
        F2 = (P_c + D_c * F1) * n;
        F3 = F1 * (p_l / ((GAMMA - 1) * r_l) - (P_c + p_l) / 2 * (1. / rho_l - 1. / r_l) + D_c * D_c / 2) + P_c * D_c;
    } else if (D_c < sn && s_r >= sn) {
        F1 = rho_r * (D_c - sn);
        F2 = (P_c + D_c * F1) * n;
        F3 = F1 * (p_r / ((GAMMA - 1) * r_r) - (P_c + p_r) / 2 * (1. / rho_r - 1. / r_r) + D_c * D_c / 2) + P_c * D_c;
    } else {
        F1 = r_r * (un_r - sn);
        F2 = (p_r + un_r * F1) * n;
        F3 = F1 * (p_r / ((GAMMA - 1) * r_r) + un_r * un_r / 2) + p_r * un_r;
    }
    dt = dx / fmax(fabs(s_l), fabs(s_r)) * CFL;

    return dt;
}
