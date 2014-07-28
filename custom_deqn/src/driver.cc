#include "driver.h"

#include <mpi.h>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "config_file.h"
#include "mesh.h"
#include "calculation.h"
#include "vtk_writer.h"
#include "data_source.h"
#include "tools.h"

Driver::Driver(const ConfigFile& config_) :  _step(0), 
                                             _config(config_) {
  MPI_Comm_rank(MPI_COMM_WORLD, &_world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &_world_size);
  parse_config();
  // Establish MPI topology
  MPI_Dims_create(_world_size, _n_dimensions, &_dim_nodes[0]);
  MPI_Cart_create(MPI_COMM_WORLD,
                  _n_dimensions,
                  &_dim_nodes[0],
                  &_dim_periods[0],
                  _mpi_reorder,
                  &_cart_comm);
  // Instantiate components of the problem
  _mesh = new Mesh(_config, _cart_comm, _dim_nodes);                            
  _calculation = new Calculation(_config, _mesh);
  _writer = new VtkWriter(_name, _mesh);
  // Set up initial conditions
  DataSource ds(_config);
  ds.populate(_mesh);
}

Driver::~Driver(){
  delete _mesh;
  delete _calculation;
  delete _writer;
}

void Driver::parse_config() {
  _debug = _config.get_or_default("debug", false);
  _name = _config.get_or_default("name", std::string("deqn_experiment"));
  std::cout << "NAME IS: " << _name << std::endl;
  _visualization_rate = _config.get_or_default("visualization_rate", 10);
  _t_start = _config.get_or_default("start_time", 0.0);                         
  _t_end = _config.get_or_default("end_time", 2.0);                             
  _del_t = _config.get_or_default("timestep", 0.02);
  if (_t_start >= _t_end) {
    throw std::logic_error("Simulation times are backwardised");
  }
  // MPI Stuff                                                                  
  _mpi_reorder = _config.get_or_default("mpi_reorder", true);                   
  _n_dimensions = _config.get_or_default("n_dimensions", 2);                 
  if (_n_dimensions != 2) {
    throw std::logic_error("Code can only handle 2D simulations at present");
  }
  _dim_nodes = _config.get_or_default("dimension_nodes", std::vector<int>());       
  _dim_periods = _config.get_or_default("dimension_periods", std::vector<int>());
  if (_dim_nodes.size() == 0) {
    _dim_nodes.resize(_n_dimensions, 0);
  } else if (_dim_nodes.size() != _n_dimensions) {
    throw std::logic_error("Dimension nodes vector is incorrectly sized");
  }
  if (_dim_periods.size() == 0) {
    _dim_periods.resize(_n_dimensions, 0);
  } else if (_dim_periods.size() != _n_dimensions) {
    throw std::logic_error("Dimension periods vector is incorrectly sized");
  }
}

void Driver::run() { 
  if (_debug) {
    std::cout << " ++ RUN COMMENCING ++ \n Initial Temperature:" 
              << Tools::mesh_sum_2d(_mesh) << std::endl;
  }
  _writer->write(0, 0.0);
  // Main application loop
  double t_now;
  for (t_now = _t_start; t_now < _t_end; t_now += _del_t) {
    _calculation->step(_del_t);
    _mesh->update_boundaries();
    _mesh->step();
    ++_step;
    if (_debug) {
      std::cout << "Step " << _step << " complete. Simulation time now " 
                << t_now << "\nTotal temp: " << Tools::mesh_sum_2d(_mesh) 
                << std::endl;
    }
    if (_step % _visualization_rate == 0 && t_now < _t_end) {
      if (_debug) { 
        std::cout << "Writing VTK file" << std::endl;
      }
      _writer->write(_step, t_now);
    }
  }
  _writer->write(_step, t_now);
  if (_debug) {
    std::cout << " ++ RUN COMPLETE ++ " << std::endl;
  }
}
