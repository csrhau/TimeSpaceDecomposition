#include "data_source.h"

#include "config_file.h"
#include "mesh.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

DataSource::DataSource(const ConfigFile& config_) : _config(config_) {
  _debug = config_.get_or_default("debug", false);
  _n_dimensions = config_.get_or_default("n_dimensions", 2);
  _dim_cells = config_.get_or_default("dimension_cells", std::vector<int>());
  _subregions = config_.get_or_default("subregions", std::vector<double>());
}
DataSource::~DataSource() {}


void DataSource::populate(Mesh * const mesh_){
  if (mesh_->n_dimensions() != 2) {
    throw std::logic_error("Incorrect dimensionality for mesh type");
  }
  double *u0 = mesh_->get_u0();
  double *u1 = mesh_->get_u1();
  int cell_count = mesh_->get_cell_count();
  // Zero initialize u1 - not strictly necessary
  std::fill(&u0[0], &u0[cell_count], 0);
  std::fill(&u1[0], &u1[cell_count], 0);
  int subregion_count = _subregions.size() / 4; // TODO This should be mesh->mdims()
  std::vector<double>::const_iterator it = _subregions.begin();
  const std::vector<int>& padded_sizes = mesh_->get_padded_sizes();
  for (int s = 0; s < subregion_count; ++s) {
    double x_min = *it++;
    double y_min = *it++;
    double x_max = *it++;
    double y_max = *it++;
    // TODO could be made faster by starting at our box.
    for (int j = 0; j < padded_sizes[1]; ++j) {
      const double y_disp = mesh_->displacement(1, j);
      for (int i = 0; i < padded_sizes[0]; ++i) {
        double x_disp = mesh_->displacement(0, i);
        if (y_disp >= y_min && y_disp < y_max
         && x_disp >= x_min && x_disp < x_max) {
          u0[i + j * padded_sizes[0]] = 10;
        }
      }
    }
  }
}
