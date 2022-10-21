#ifndef GODUNOV_SOLVER_SYSTEM_HELPER_HPP
#define GODUNOV_SOLVER_SYSTEM_HELPER_HPP

#include <vector>
#include <string>

auto mk_dir(const std::vector<double>& parameters) -> std::string;

void plots(const std::string& file_name);

#endif //GODUNOV_SOLVER_SYSTEM_HELPER_HPP
