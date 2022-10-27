#include <iostream>
#include "../gas_parameters.hpp"

std::ostream &operator<<(std::ostream &os, const gas_parameters &gas){
  os << std::scientific << gas.r << "\t" << gas.u << "\t" << gas.p << "\t";
  return os;
}