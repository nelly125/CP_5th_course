#include <stdexcept>
#include <iostream>
#include <sys/stat.h>
#include <sstream>
#include <iomanip>
#include "system_helper.hpp"

#include <filesystem>

auto mk_dir( uint32_t N, double time ) -> std::string {

  std::cout << std::filesystem::current_path() << std::endl;

  std::string directory = "./results/piston_";
  uint32_t check;

  std::ostringstream streamObj3;
  streamObj3 << std::fixed;
  streamObj3 << std::setprecision(2);
  streamObj3 << time;

  std::string temp_string = std::to_string(N);
  directory += temp_string + "_";
  directory += streamObj3.str();
/*  temp_string = std::to_string(SIGMA);
  directory += temp_string + "_";
  temp_string = std::to_string(AMPLITUDE);
  directory += temp_string + "_";
  temp_string = std::to_string(OMEGA);
  directory += temp_string + "_";*/

  std::cout << directory << std::endl;

  check = mkdir(directory.c_str(), 0777);
  if (!check)
    std::cout << "Directory created\n";
  else {
    std::cout << "Unable to create directory\n";
  }
  directory += "/";

  std::string plot_dir = directory + "plots/";

  check = mkdir(plot_dir.c_str(), 0777);
  if (!check)
    std::cout << "Directory created\n";
  else {
    std::cout << "Unable to create directory\n";
  }

  return directory;
}

void plots( const std::string &directory, uint32_t N ) {
  std::string temp_string;
  std::string python_script;

  temp_string = std::to_string(N);
  python_script = "python3 ./src/plots_scripts/animated_plots.py ";
  python_script += temp_string + " ";
  python_script += directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }

  temp_string = std::to_string(N);
  python_script = "python3 ./src/plots_scripts/trajectories_plots.py ";
  python_script += temp_string + " ";
  python_script += directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }

  temp_string = std::to_string(N);
  python_script = "python3 ./src/plots_scripts/delta_energy_plots.py ";
  python_script += temp_string + " ";
  python_script += directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }

  temp_string = std::to_string(N);
  python_script = "python3 ./src/plots_scripts/energy_plots.py ";
  python_script += temp_string + " ";
  python_script += directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }
}