#ifndef DRIVER_H
#define DRIVER_H

#include <mpi.h>
#include <string>
#include <vector>


class ConfigFile;
class VtkWriter;
class Mesh;
class Calculation;

class Driver {
 public:
  Driver(const ConfigFile& config_);
  ~Driver();
  void run();

 private:
  void parse_config();
  bool _debug;
  int _step;
  int _visualization_rate;
  double _t_start;
  double _t_end;
  double _del_t;
  size_t _n_dimensions;
  std::string _name;
  const ConfigFile& _config;
  Mesh * _mesh;
  Calculation * _calculation;
  VtkWriter * _writer;
  // MPI members
  int _world_rank;
  int _world_size;
  bool _mpi_reorder;
  std::vector<int> _dim_nodes;
  std::vector<int> _dim_periods;
  MPI_Comm _cart_comm;
};
#endif
