#include "data_source.h"

#include "config_file.h"
#include "mesh.h"

#include <vector>
#include <algorithm>
#include <stdexcept>

DataSource::DataSource(const ConfigFile& config_) : _config(config_) {
  _debug = config_.get_or_default("debug", false);
  _subregions = config_.get_or_default("subregions", std::vector<double>());
}
DataSource::~DataSource() {}

void DataSource::populate(Mesh * const mesh_){
  double *u0 = mesh_->get_u0();
  double *u1 = mesh_->get_u1();
  int augmented_cell_count = mesh_->get_node_augmented_cell_count();
  // Zero initialize u1 - not strictly necessary
  std::fill(&u0[0], &u0[augmented_cell_count], 0);
  std::fill(&u1[0], &u1[augmented_cell_count], 0);
  int subregion_count = _subregions.size() / 4; // x, y * 2 dims
  std::vector<double>::const_iterator it = _subregions.begin();
  for (int s = 0; s < subregion_count; ++s) {
    double x_min = *it++;
    double y_min = *it++;
    double x_max = *it++;
    double y_max = *it++;
    for (int i = 0; i < mesh_->get_node_augmented_row_count(); ++i) {
      double x_coord = mesh_->get_y_coord(i);
      if (x_min <= x_coord && x_coord < x_max) {
        for (int j = 0; j < mesh_->get_node_augmented_col_count(); ++j) {
          double y_coord = mesh_->get_x_coord(j); 
          if (y_min <= y_coord && y_coord < y_max) {
            u0[i * mesh_->get_node_augmented_col_count() + j] = 10;
          }
        }
      }
    }
  }
}
