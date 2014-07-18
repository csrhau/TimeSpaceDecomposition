#include "experiment.h"
#include <iostream>
#include "config_file.h"

Experiment::Experiment(const ConfigFile config_) : _config(config_) {
  bool debug;
  std::cout << "Experiment starting" << std::endl;

}

Experiment::~Experiment() {}


void Experiment::run() {
  
}
