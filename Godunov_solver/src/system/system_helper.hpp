#ifndef GODUNOV_SOLVER_SYSTEM_HELPER_HPP
#define GODUNOV_SOLVER_SYSTEM_HELPER_HPP

#include <vector>
#include <string>

//auto mk_dir(const std::vector<double>& parameters) -> std::string;

auto mk_dir(int N, double time) -> std::string;

void plots(const std::string &directory, int N);

#endif //GODUNOV_SOLVER_SYSTEM_HELPER_HPP
