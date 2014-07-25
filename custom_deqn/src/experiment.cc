#include "experiment.h"
#include <iostream>
#include <vector>
#include "config_file.h"
#include "data_source.h"
#include "calculation.h"

#include "mesh.h"
#include "vtk_writer.h"
#include "tools.h"

namespace {

} // File-scoped

Experiment::Experiment(const ConfigFile * const config_) : _config(config_),
                                                           _step(0) {
  _debug = _config->get_or_default("debug", false);
  _name = _config->get_or_default("name", _config->get_filename());
  std::cout << "Using " << _name << " as experiment name" << std::endl;
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
  VtkWriter writer(_name, _mesh);
  std::cout << "experiment running!!" << std::endl;
  if (_debug) {
    std::cout << "+++++++++++++++++++++++++++++\n"
                 "+ RUN COMMENCING - STAND BY +\n"
                 "+++++++++++++++++++++++++++++" << std::endl;
    std::cout << "Total temp: " << Tools::mesh_sum_2d(_mesh) << "\n";
  }
  // Initial dump
  writer.write(0, 0.0);
  // Main application Loop
  double t_now;
  for (t_now = _t_start; t_now < _t_end; t_now += _del_t) {
    _calculation->step(_del_t);
    _mesh->update_boundaries();
    _mesh->step();
    ++_step;
    if (_step % _visualization_rate == 0 && t_now < _t_end) {
      if (_debug) {
        std::cout << "Total temp: " << Tools::mesh_sum_2d(_mesh) << "\n";
      }
      writer.write(_step, t_now);
    }
    if (_debug) {
      std::cout << "Step " << _step << " complete. Total temp: " 
                << Tools::mesh_sum_2d(_mesh) << "\n";
    }
  }
  // Final dump
  writer.write(_step, t_now);
  if (_debug) {
    std::cout << "+++++++++++++++++++++++++++++\n"
                 "+ RUN FINISHED SUCCESSFULLY +\n"
                 "+++++++++++++++++++++++++++++" << std::endl;
  }
}
