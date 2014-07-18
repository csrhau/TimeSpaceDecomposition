#include <iostream>
#include <string>

#include "config_file.h"
#include "experiment.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: deqn <filename>" << std::endl;
    return 1;
  }
  std::string filename = std::string(argv[1]);
  const ConfigFile config(filename);
  Experiment experiment(config);
  experiment.run();
  std::cout << "Experiment complete!" << std::endl;
  return 0;
}
