#include <stdexcept>
#include <iostream>
#include <sys/stat.h>
#include <iomanip>
#include "system_helper.hpp"

#include <filesystem>

auto mk_dir( uint32_t N, double time, double amplitude, double omega, double sigma, double diaph ) -> std::string {

  std::cout << std::filesystem::current_path() << std::endl;

  std::string directory = "./results/piston__";
  uint32_t check;

  std::ostringstream streamObj3;
  streamObj3 << std::fixed;
  streamObj3 << std::setprecision(2);
  streamObj3 << time << "__" << amplitude << "__" << omega << "__" << sigma << "__" << diaph;

  std::string temp_string = std::to_string(N);
  directory += temp_string + "__";
  directory += streamObj3.str();


/*  temp_string = std::to_string(amplitude);
  streamObj3 << amplitude;
  directory += temp_string + "_";
  temp_string = std::to_string(omega);
  directory += temp_string + "_";
  temp_string = std::to_string(sigma);
  directory += temp_string + "_";
  temp_string = std::to_string(diaph);
  directory += temp_string + "_";*/

  std::cout << directory << std::endl;

  check = mkdir(directory.c_str(), 0777);
  if (!check)
    std::cout << "Directory created\n";
  else {
    std::cout << "Unable to create directory\n";
  }
  directory += "/";

  std::string data_dir = directory + "data/";
  check = mkdir(data_dir.c_str(), 0777);
  if (!check)
    std::cout << "Directory created\n";
  else {
    std::cout << "Unable to create directory\n";
  }

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
  python_script = "python3 ./src/python_scripts/run_all_scripts.py ";
  python_script += temp_string + " ";
  python_script += directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }

/*  temp_string = std::to_string(N);
  python_script = "python3 ./src/python_scripts/animated_plots.py ";
  python_script += temp_string + " ";
  python_script += directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }

  temp_string = std::to_string(N);
  python_script = "python3 ./src/python_scripts/trajectories_plots.py ";
  python_script += temp_string + " ";
  python_script += directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }

  temp_string = std::to_string(N);
  python_script = "python3 ./src/python_scripts/delta_energy_plots.py ";
  python_script += temp_string + " ";
  python_script += directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }

  temp_string = std::to_string(N);
  python_script = "python3 ./src/python_scripts/energy_plots.py ";
  python_script += temp_string + " ";
  python_script += directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }*/
}