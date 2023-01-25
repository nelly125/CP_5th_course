#ifndef GODUNOV_SOLVER_SYSTEM_HELPER_HPP
#define GODUNOV_SOLVER_SYSTEM_HELPER_HPP

#include <vector>
#include <string>

//auto mk_dir(const std::vector<double>& parameters) -> std::string;

auto mk_dir( uint32_t N, double time, double amplitude, double omega, double sigma, double diaph ) -> std::string;

auto mk_dir( uint32_t N, double time) -> std::string;

void plots( const std::string &directory, uint32_t N );

void spherical_plots( const std::string &directory);

#endif //GODUNOV_SOLVER_SYSTEM_HELPER_HPP
