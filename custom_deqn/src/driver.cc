#include "driver.h"

#include <iostream>
#include <vector>
#include <stdexcept>

#include "config_file.h"
#include "data_source.h"
#include "calculation.h"

#include "mesh.h"
#include "vtk_writer.h"
#include "tools.h"

Driver::Driver(const ConfigFile& config_) : _config(config_),
                                            _step(0) {
  parse_config();
  init();
}

Driver::~Driver() {
  delete _calculation;
  delete _mesh;
}

void Driver::parse_config(){
  _debug = _config.get_or_default("debug", false);
  _name = _config.get_or_default("name", _config.get_filename());
  std::cout << "Using " << _name << " as experiment name" << std::endl;
  _visualization_rate = _config.get_or_default("visualization_rate", 5);
  _t_start = _config.get_or_default("start_time", 0.0);
  _t_end = _config.get_or_default("end_time", 2.0);
  _del_t = _config.get_or_default("timestep", 0.02); // TODO...variable dt?
  // MPI Stuff
  _mpi_reorder = _config.get_or_default("mpi_reorder", true);
  _n_dimensions = _config.get_or_default("n_dimensions", 2);  
  _dim_nodes= _config.get_or_default("dimension_nodes", std::vector<int>());
  if (_dim_nodes.size() == 0) {
    for (size_t i = 0; i < _n_dimensions; ++i) {
      _dim_nodes.push_back(0);
    }
  }
  _dim_periods = _config.get_or_default("dimension_periods", std::vector<int>());
  if (_dim_periods.size() == 0) {
    for (size_t i = 0; i < _n_dimensions; ++i) {
      _dim_periods.push_back(0);
    }
  }
  if (_dim_periods.size() != _n_dimensions || _dim_nodes.size() != _n_dimensions) {
    throw std::logic_error("Period/dimension information does mismatch");
  }
}

void Driver::init() {
  MPI_Comm_rank(MPI_COMM_WORLD, &_world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &_world_size);
  MPI_Dims_create(_world_size, _n_dimensions, &_dim_nodes[0]);
  std::cout << "Debug: " << _debug << " world rank " << _world_rank << std::endl;
  if (_debug &&(_world_rank == 0)) {
    std::cout << "Creating a ";
    size_t d;
    for (d = 0; d < _n_dimensions-1; ++d) {
      std::cout << _dim_nodes[d] << " x ";
    }
    std::cout << _dim_nodes[d] << " space" << std::endl;
  }

  MPI_Cart_create(MPI_COMM_WORLD,
                  _n_dimensions,
                  &_dim_nodes[0],
                  &_dim_periods[0],
                  _mpi_reorder,
                  &_cart_comm);

  _mesh = new Mesh(_config, _cart_comm, _dim_nodes);
  _calculation = new Calculation(_config, _mesh);
  DataSource ds(_config);
  ds.populate(_mesh);
}

void Driver::run() {
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
