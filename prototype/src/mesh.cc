#include "mesh.h"

#include <vector>
#include "config_file.h"

Mesh::Mesh(const ConfigFile& config_) : _config(config_) {
  // Calculate our simulation domain
  // core space excludes ghost cells and boundary padding
  std::vector<int> core_dimensions = _config.get_or_default("logical_dimensions",
                                                             std::vector<int>());
  _world_core_row_count = core_dimensions.at(0);
  _world_core_col_count = core_dimensions.at(1);
  std::vector<double> physical_dimensions = _config.get_or_default("physical_dimensions",
                                                             std::vector<double>());
  _world_height = physical_dimensions.at(0);
  _world_width = physical_dimensions.at(1);
}

double Mesh::get_world_core_row_count() const {
 return _world_core_row_count;
}
double Mesh::get_world_core_col_count() const {
 return _world_core_col_count;
}
double Mesh::get_world_height() const {
 return _world_height;
}
double Mesh::get_world_width() const {
 return _world_width;
}
double Mesh::get_del_y() const { 
  return _world_height / _world_core_row_count; 
}
double Mesh::get_del_x() const { 
  return _world_width / _world_core_col_count; 
}
