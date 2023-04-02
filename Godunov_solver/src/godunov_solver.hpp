#ifndef GODUNOV_SOLVER_GODUNOV_SOLVER_HPP
#define GODUNOV_SOLVER_GODUNOV_SOLVER_HPP

#include <fstream>
#include <vector>
#include <functional>
#include "../gas_parameters.hpp"

enum boundaries { piston, soft, hard_wall };

class solver {
public:
  solver() = default;

  ~solver() = default;
/*
    void solve_system(double x_0, double x_n, double x_c, double time, uint32_t N,
                      double r, double u, double p);*/

  static std::string solve_system( double x_0,
                                   double x_n,
                                   double x_c,
                                   double time,
                                   uint32_t N,
                                   double r_l,
                                   double u_l,
                                   double p_l,
                                   double r_r,
                                   double u_r,
                                   double p_r,
                                   double _amplitude,
                                   double _omega,
                                   const std::function<double( double x_0,
                                                               double left_0,
                                                               double amplitude,
                                                               double omega,
                                                               double time )> &P_func,
                                   const std::function<double( double x_0,
                                                               double left_0,
                                                               double amplitude,
                                                               double omega,
                                                               double time )> &right_piston_func,
                                   boundaries left_boundary,
                                   boundaries right_boundary,
                                   bool hllc_flag );

  static std::string solve_system_boundaries( double x_0,
                                              double x_n,
                                              double time,
                                              uint32_t N,
                                              double r_,
                                              double u_,
                                              double p_,
                                              double _amplitude,
                                              double _omega,
                                              const std::function<double( double,
                                                                   double,
                                                                   double,
                                                                   double,
                                                                   double )> &left_piston_func,
                                              const std::function<double( double,
                                                                   double,
                                                                   double,
                                                                   double,
                                                                   double )> &right_piston_func,
                                              boundaries left_boundary,
                                              boundaries right_boundary,
                                              bool hllc_flag );
private:

  static double find_s_cell( gas_parameters left, gas_parameters right );
  static double find_x( uint32_t j, double left, double diaph, double right, uint32_t N, uint32_t i_contact );
  static double compute_dt( const std::vector<gas_parameters> &,
                            uint32_t N,
                            uint32_t i_contact,
                            double dx_left_0,
                            double dx_right_0 );

  static void output_parameters_to_file( std::ostream &out,
                                         const std::vector<gas_parameters> &gas,
                                         double time,
                                         double left,
                                         double diaph,
                                         double right,
                                         uint32_t i_contact,
                                         uint32_t N );

  static void output_parameters_to_file( std::ostream &out,
                                         const std::vector<gas_parameters> &gas,
                                         double time,
                                         double left,
                                         double right,
                                         uint32_t N );

  static void trajectory_to_file( std::ostream &out,
                                  double time,
                                  double left,
                                  double diaph,
                                  double right );

  static double total_energy( const std::vector<gas_parameters> &gas,
                              uint32_t i_contact,
                              double dx_left,
                              double dx_right );
  static double total_kinetic_energy( const std::vector<gas_parameters> &gas );
  static double total_internal_energy( const std::vector<gas_parameters> &gas );
  static void output_cell_params_to_file( std::ostream &out,
                                          const std::vector<gas_parameters> &gas,
                                          double time,
                                          double left,
                                          double diaph,
                                          double right,
                                          uint32_t i_contact,
                                          uint32_t N );
  static double total_energy( const std::vector<gas_parameters> &gas, double dx );
  static void trajectory_to_file( std::ostream &out, double time, double left, double right );
  static double compute_dt( const std::vector<gas_parameters> &gas, uint32_t N, double dx );
};

#endif //GODUNOV_SOLVER_GODUNOV_SOLVER_HPP
