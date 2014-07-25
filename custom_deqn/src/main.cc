#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include "config_file.h"
#include "experiment.h"

bool check_config_(const ConfigFile* config_) {
  bool valid = true;
  if (config_->get_or_default("debug", false)) {
    config_->print_config();
  }

  std::size_t n_dims = config_->get_or_default("n_dimensions", 0);
  if (n_dims <= 0) {
    std::cerr << "Error! n_dimensions has to be greater than zero, currently:"
              << n_dims << std::endl;
    valid = false;
  }


  std::vector<int> dim_cells = config_->get_or_default("dimension_cells",
                                                      std::vector<int>());
  if (dim_cells.size() == 0) {
    std::cerr << "Error! Config is missing dimension_cells key:\n"
                 "\t e.g. dimension_cells 500 400" << std::endl;
    valid = false;
  } else if (dim_cells.size() != n_dims) {
    std::cerr << "Error! the length of the dimension_cells array does not "
                 "match the n_dimensions parameter";
    valid = false;
  }


  std::vector<double> min_coords =
    config_->get_or_default("min_coordinates", std::vector<double>());
  std::vector<double> max_coords =
    config_->get_or_default("max_coordinates", std::vector<double>());

  if (min_coords.size() == 0 || max_coords.size() == 0) {
    std::cerr << "Have to specify the coordinates of the minimum and maximum "
                 "points in the space" << std::endl;
    valid = false;
  } else if (min_coords.size() != n_dims || max_coords.size() != n_dims) {
    std::cerr << "Min/Max coordinates have to be specified for the correct "
                 "number of dimensions!" << std::endl;
    valid = false;
  }








  std::vector<double> regions = config_->get_or_default("subregions",
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
  const ConfigFile config_(filename);


  if (!check_config_(&config_)) {
    std::cout << "Bad config! exiting" << std::endl;
    return 1;
  }


  try {
    std::cout << "Construction of experiment" << std::endl;
    Experiment experiment(&config_);
    std::cout << "Construction of experiment done " << std::endl;
    experiment.run();
    std::cout << "Experiment complete!" << std::endl;
  } catch (std::logic_error& ex) {
    std::cerr << "Exception thrown!! " << ex.what();
    return 1;
  }
  return 0;
}
