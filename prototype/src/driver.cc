#include "driver.h"

#include <mpi.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "tools-inl.h"
#include "config_file.h"
#include "vtk_writer.h"
#include "data_source.h"
#include "mesh.h"
#include "static_mesh.h"
#include "static_blocking_mesh.h"
#include "calculation.h"

Driver::Driver(const ConfigFile& config_) : _config(config_) {
  // Read configuration
  _debug = _config.get_or_default("debug", false);
  _name = _config.get_or_default("name", std::string("prototype"));
  _visualization_rate = _config.get_or_default("visualization_rate", 1);
  _t_start = _config.get_or_default("start_time", 0.0);
  _t_end = _config.get_or_default("end_time", 2.0);
  _del_t = _config.get_or_default("timestep", 0.02);
  _mpi_reorder = _config.get_or_default("mpi_reorder", true);
  _dim_nodes = _config.get_or_default("dim_nodes", std::vector<int>());       
  _dim_periods = _config.get_or_default("dim_periods", std::vector<int>());      
  // Resize to be 2 elements (1 per dimension), 0 pad
  _dim_nodes.resize(2, 0); 
  _dim_periods.resize(2, 0);
  // Establish MPI topology
  MPI_Comm_size(MPI_COMM_WORLD, &_world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &_world_rank);
  MPI_Dims_create(_world_size, 2/*2D*/, &_dim_nodes[0]); 
  MPI_Cart_create(MPI_COMM_WORLD,
                  2, // 2D
                  &_dim_nodes[0],
                  &_dim_periods[0],
                  _mpi_reorder,
                  &_cart_comm);
  // Initialize Mesh (and calculation)
  _mesh_type = _config.get_or_default("mesh_type", std::string("static"));
  if (_mesh_type == "static") {
    _mesh = new StaticMesh(_config, _cart_comm, _dim_nodes);
  } else if (_mesh_type == "static_blocking") {
    _mesh = new StaticBlockingMesh(_config, _cart_comm, _dim_nodes);
  } else {
    std::stringstream ss;
    ss << "Unknown mesh type: " << _mesh_type << std::endl;
    throw std::logic_error(ss.str());
  }
  std::stringstream ss;
  ss << _name << "_" << _mesh_type;
  _outfile_tag = ss.str();

  _calculation = new Calculation(_config, _mesh);
  // Datasource initialize
  DataSource ds(_config);
  ds.populate(_mesh);
}

Driver::~Driver(){
  delete _mesh;
  delete _calculation;
}

void Driver::run() {
  VtkWriter writer(_outfile_tag, _mesh, _world_rank, _world_size);
  double wall_start, wall_stop;
  double cpu_start, cpu_stop;
  if (_debug) {
    std::cout << " ++ RUN BEGINNING ++ " << std::endl;
  }
  timers(wall_start, cpu_start); // start timing
  int step = 0;
  double t_now = _t_start;
  while (t_now < _t_end) { // doublecompare
    if (step % _visualization_rate == 0) {
      writer.write(step, t_now);
      if (_debug) {
        double temp = local_temp();
        double global_temp = 0;
        MPI_Reduce(&temp, &global_temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (_world_rank == 0) {
        std::cout << " Outputting vtk file for step " << step << ",\n\ttnow = "
                  << t_now << ",\n\tvis rate:" << _visualization_rate 
                  << "\n\ttotal temp:" << global_temp <<  std::endl;
        }
      }
    }
    _calculation->step(_del_t);
    _mesh->advance();
    ++step;
    t_now += _del_t;
  }
  writer.write(step, t_now);
  timers(wall_stop, cpu_stop); // stop timing
  if (_debug) {
    std::cout << " ++ RUN FINISHING ++ " << std::endl;
  }
  // Output timing information
  MPI_Barrier(MPI_COMM_WORLD);
  if (_world_rank == 0) {
    std::cout << "Timings: wallclock:" << (wall_stop - wall_start) << "s\n"
                 "         cpu clock:" << (cpu_stop - cpu_start) << std::endl;
  }
}

double Driver::local_temp() const {
  double total = 0;
  int x_span = _mesh->get_node_augmented_col_count();
  double *u0 = _mesh->get_u0();
  for (int i = 1; i < _mesh->get_node_core_row_count() + 1; ++i) {
    for (int j = 1; j < _mesh->get_node_core_col_count() + 1; ++j) {
      total += u0[i * x_span + j];
    }
  }
  return total;
}
