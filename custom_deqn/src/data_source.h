#ifndef DATA_SOURCE_H
#define DATA_SOURCE_H

#include <vector>

class ConfigFile;
class Mesh;
class DataSource {
 public: 
  DataSource(const ConfigFile * const config_);
  ~DataSource();
  void populate(Mesh * const mesh);

 private:
  const ConfigFile * const  _config;
  bool _debug;
  int _n_dimensions;
  std::vector<int> _dim_sizes;
  std::vector<int> _subregions;

    
  int _nx;
  int _ny;
  void init();
};

#endif
