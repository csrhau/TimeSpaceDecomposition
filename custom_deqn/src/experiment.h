#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "calculation.h"

class ConfigFile;
class Mesh;
class Experiment {
 public:
  Experiment(const ConfigFile * const config_);
  ~Experiment();
  void run();

 private:
  const ConfigFile * const  _config;
  int _step;
  Calculation* _calculation; // todo const ptrs
  Mesh* _mesh;
  bool _debug;
  int _visualization_rate;
  double _t_start;
  double _t_end;
  double _del_t;

}; 

#endif
