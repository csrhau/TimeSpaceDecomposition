#include "experiment.h"
#include <sys/time.h>
#include <iostream>
#include "config_file.h"
#include "data_source.h"
#include "calculation.h"
#include "mesh.h"

Experiment::Experiment(const ConfigFile * const config_) : _config(config_),
                                                           _step(0) {
  std::cout << "Experiment constructor " << __LINE__ << std::endl;
  _debug = _config->get_or_default("debug", false);
  _visualization_rate = _config->get_or_default("visualization_rate", 5);
  _t_start = _config->get_or_default("start_time", 0.0);
  _t_end = _config->get_or_default("end_time", 2.0);
  _del_t = _config->get_or_default("timestep", 0.02); // TODO...variable dt?
  _mesh = new Mesh(_config);
  _calculation = new Calculation(_config, _mesh);
  DataSource ds(_config);
  ds.populate(_mesh);
}

Experiment::~Experiment() {
  delete _calculation;
  delete _mesh;
}

void Experiment::run() {
  std::cout << "experiment running!!" << std::endl;
  struct timeval start_time, end_time;
  if (_debug) {
    std::cout << "+++++++++++++++++++++++++++++\n"
                 "+ RUN COMMENCING - STAND BY +\n"
                 "+++++++++++++++++++++++++++++" << std::endl;
    gettimeofday(&start_time, NULL);
  }
  // Write out initial problem
  std::cout << "TODO: Write out initial step " << _step << std::endl;
  // Main application LOOP
  for (double t_now = _t_start; t_now < _t_end; t_now += _del_t) {
    // Perform calculation and updates
    _calculation->step();
    _mesh->step();
    ++_step;
    // Produce required outputs
    if (_step % _visualization_rate == 0 && t_now < _t_end) {
      std::cout << "TODO: Write out intermediate step " << _step << std::endl;
    }
    if (_debug) {
      std::cout << "Step " << _step << " complete" << std::endl;
    }
  }
  std::cout << "TODO: Write out final step " << _step << std::endl;


  if (_debug) {
    std::cout << "+++++++++++++++++++++++++++++\n"
                 "+ RUN FINISHED SUCCESSFULLY +\n"
                 "+++++++++++++++++++++++++++++" << std::endl;
    gettimeofday(&end_time, NULL);
    // TODO finish timings

  }
}


