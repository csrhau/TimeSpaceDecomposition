#ifndef DATA_SOURCE_H
#define DATA_SOURCE_H

#include <vector>

class ConfigFile;
class Mesh;
class DataSource {
 public: 
  DataSource(const ConfigFile * const config_);
  ~DataSource();
  void populate(Mesh * const mesh_);

 private:
  const ConfigFile * const  _config;
  bool _debug;
  int _n_dimensions;
  std::vector<int> _dim_cells;
  std::vector<double> _subregions;

    
  int _nx;
  int _ny;
  void init();
};

#endif
