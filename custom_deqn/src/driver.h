#ifndef DRIVER_H
#define DRIVER_H

#include "calculation.h"
#include <string>

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
  std::string _name;

};

#endif
