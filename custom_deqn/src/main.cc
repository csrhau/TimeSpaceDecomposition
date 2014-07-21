#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include "config_file.h"
#include "experiment.h"

bool check_config(const ConfigFile* config) {
  bool valid = true;
  if (config->get_or_default("debug", false)) {
    config->print_config();
  }
  if (!config->get_or_default("dimension_sizes", NULL)) {
    std::cerr << "Error! Config is missing dimension_sizes key:\n"
                 "\t e.g. dimension_sizes 500 400" << std::endl;
    valid = false;
  }

  std::vector<double> regions = config->get_or_default("subregions",
                                                    std::vector<double>());
  if (regions.size() % 4 != 0) {
    std::cerr << "Error! Subregions have to be square!\n"
                 "\t e.g. subregions  10 10 30 30" << std::endl;
    valid = false;
  }
  return valid;
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: deqn <filename>" << std::endl;
    return 1;
  }
  std::string filename = std::string(argv[1]);
  const ConfigFile config(filename);


  if (!check_config(&config)) {
    std::cout << "Bad config; exiting" << std::endl;
    return 1;
  }


  try {
    std::cout << "Construction of experiment" << std::endl;
    Experiment experiment(&config);
    std::cout << "Construction of experiment done " << std::endl;
    experiment.run();
    std::cout << "Experiment complete!" << std::endl;
  } catch (std::logic_error& ex) {
    std::cerr << "Exception thrown!! " << ex.what();
    return 1;
  }
  return 0;
}
