#ifndef DRIVER_H
#define DRIVER_H

#include <mpi.h>
#include <vector>
#include <string>

class ConfigFile;
class Mesh;
class Calculation;

class Driver {
 public:
  Driver(const ConfigFile& config_);
  ~Driver();
  void run();

 private:
  double local_temp() const;
  bool _debug;
  std::string _name;
  int _visualization_rate;
  double _t_start;
  double _t_end;
  double _del_t;
  const ConfigFile& _config;
  Mesh * _mesh;
  Calculation * _calculation;
  // MPI members
  std::vector<int> _dim_nodes;
  std::vector<int> _dim_periods;
  MPI_Comm _cart_comm;
  bool _mpi_reorder;
  int _world_size;
  int _world_rank;
  int _cart_size;
  int _cart_rank;


};
#endif
