#include "experiment.h"
#include <iostream>
#include <vector>
#include "config_file.h"
#include "data_source.h"
#include "calculation.h"
#include "mesh.h"

namespace {
// TODO move into common header file
#define POLY2(i, j, ispan) ((i) + ((j)  * (ispan)))
void print_mesh(Mesh* _mesh) {
  double* space = _mesh->get_u0();
  const std::vector<int>& padded_sizes = _mesh->get_padded_sizes();
  for (int j = 0; j < padded_sizes[1]; ++j) {
    for (int i = 0; i < padded_sizes[0]; ++i) {
      std::cout << "\t" <<  space[i + j * padded_sizes[0]]<<",";
    }
    std::cout << std::endl;
  }
}

static void print_total_temp(Mesh* _mesh) {
  double *u0 = _mesh->get_u0();
  int x_min = _mesh->get_from_index(0);
  int y_min = _mesh->get_from_index(1);
  int x_max = _mesh->get_to_index(0);
  int y_max = _mesh->get_to_index(1);
  int x_span = _mesh->get_dimension_span(0);
  double temp = 0;
  for (int j = y_min; j < y_max; ++j) { //TODO deal with boundary
    for (int i = x_min; i < x_max; ++i) {
      int center = POLY2(i, j, x_span);
      temp += u0[center];
    }
  }
  std::cout << "Total temperature: " << temp << std::endl;
}
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
  std::cout << "experiment running!!" << std::endl;
  if (_debug) {
    std::cout << "+++++++++++++++++++++++++++++\n"
                 "+ RUN COMMENCING - STAND BY +\n"
                 "+++++++++++++++++++++++++++++" << std::endl;
    print_mesh(_mesh);
  }
  // Main application Loop
  for (double t_now = _t_start; t_now < _t_end; t_now += _del_t) {
    _calculation->step(_del_t);
    _mesh->update_boundaries();
    _mesh->step();
    ++_step;
    if (_step % _visualization_rate == 0 && t_now < _t_end) {
      if (_debug) print_mesh(_mesh);
      std::cout << "TODO: Output visualization file" << std::endl;
    }
    if (_debug) {
      std::cout << "Step " << _step << " complete" << std::endl;
      print_total_temp(_mesh);

    }
  }
  if (_debug) {
    std::cout << "+++++++++++++++++++++++++++++\n"
                 "+ RUN FINISHED SUCCESSFULLY +\n"
                 "+++++++++++++++++++++++++++++" << std::endl;
  }
}
