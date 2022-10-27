#ifndef GODUNOV_SOLVER_GODUNOV_SOLVER_HPP
#define GODUNOV_SOLVER_GODUNOV_SOLVER_HPP

#include <fstream>
#include <vector>
#include "../gas_parameters.hpp"
class solver {
public:
    solver() = default;

    ~solver() = default;
/*
    void solve_system(double x_0, double x_n, double x_c, double time, int N,
                      double r, double u, double p);*/

    void solve_system(double x_0, double x_n, double x_c, double time, int N,
                      double r_l, double u_l, double p_l,
                      double r_r, double u_r, double p_r);

private:

  void output_parameters_to_file(std::ostream out, const std::string& dir);
  static double find_s_cell(gas_parameters left, gas_parameters right);
  static double find_x(int j, double left, double diaph, double right, int N, int i_contact);
  static void get_values(int k, int N, bool minmod_flag, gas_parameters gas, double &r, double &u, double &p);
  static double compute_dt(const std::vector<gas_parameters>&,
                    int N,
                    int i_contact,
                    double dx_left_0,
                    double dx_right_0);
  static void output_parameters_to_file(std::ostream& out, const std::vector<gas_parameters>& gas, double time);
  void output_parameters_to_file(std::ostream &out,
                                 const std::vector<gas_parameters> gas,
                                 double time,
                                 double left,
                                 double diaph,
                                 double right,
                                 int i_contact,
                                 int N);
};

#endif //GODUNOV_SOLVER_GODUNOV_SOLVER_HPP
