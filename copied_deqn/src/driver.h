#ifndef DRIVER_H
#define DRIVER_H

#include "calculation.h"

#include <vector>
#include <string>
#include <mpi.h>


class ConfigFile;
class Mesh;
class Driver {
 public:
  Driver(const ConfigFile& config_);
  ~Driver();
  void run();

 private:
  const ConfigFile& _config;
  int _step;
  Calculation *_calculation;
  Mesh *_mesh;
  bool _debug;
  int _visualization_rate;
  double _t_start;
  double _t_end;
  double _del_t;
  size_t _n_dimensions;
  std::string _name;
  // MPI Things
  bool _mpi_reorder;
  int _world_rank;
  int _world_size;
  std::vector<int> _dim_nodes;
  std::vector<int> _dim_periods;
  MPI_Comm _cart_comm;
  void parse_config();
  void init();

};

#endif
