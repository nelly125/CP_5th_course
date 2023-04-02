#include <iostream>
#include "../gas_parameters.hpp"

std::ostream &operator<<( std::ostream &os, const gas_parameters &gas ) {
  os << gas.r << "\t" << gas.u << "\t" << gas.p << "\t";
  return os;
}

flows gas_to_flows(gas_parameters &gas) {
  double e = gas.p / ((GAMMA - 1) * gas.r);
  flows flow{};
  flow.M = gas.r * gas.u;
  flow.J = gas.p + gas.r * gas.u * gas.u;
  flow.E = gas.r * gas.u * (e + (0.5 * gas.u * gas.u)) + gas.p * gas.u;
  return flow;
}


cons_flows gas_to_cons_flow(gas_parameters &gas) {
  cons_flows c_flow{};
  double e = gas.p / ((GAMMA - 1) * gas.r);
  c_flow.m = gas.r;
  c_flow.j = gas.r * gas.u;
  c_flow.e = gas.r * (e + (0.5 * gas.u * gas.u));
  return c_flow;
}

gas_parameters c_flow_to_gas(cons_flows &c_flow) {
  gas_parameters gas{};
  gas.r = c_flow.m;
  gas.u = c_flow.j / c_flow.m;
  gas.p = (GAMMA - 1) * (c_flow.e - ((0.5 * c_flow.j * c_flow.j) / c_flow.m));
  return gas;
}


flows flows::operator-(flows tmp) {
  flows res{};
  res.M = M - tmp.M;
  res.J = J - tmp.J;
  res.E = E - tmp.E;
  return res;
}