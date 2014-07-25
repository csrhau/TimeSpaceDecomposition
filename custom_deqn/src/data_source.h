#ifndef DATA_SOURCE_H
#define DATA_SOURCE_H

#include <vector>

class ConfigFile;
class Mesh;
class DataSource {
 public:
  DataSource(const ConfigFile& config_);
  ~DataSource();
  void populate(Mesh * const mesh_);

 private:
  const ConfigFile& _config;
  bool _debug;
  int _n_dimensions;
  std::vector<int> _dim_cells;
  std::vector<double> _subregions;

  int _nx;
  int _ny;
};

#endif
