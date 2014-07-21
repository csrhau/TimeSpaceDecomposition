#include "data_source.h"

#include "config_file.h"
#include "mesh.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace {
struct RegionBoundaries {
  double xmin;
  double ymin;
  double xmax;
  double ymax;
};

}

DataSource::DataSource(const ConfigFile * const config_) : _config(config_) {
  _debug = config_->get_or_default("debug", false);
  _n_dimensions = config_->get_or_default("n_dimensions", 2);
  _dim_sizes = config_->get_or_default("dimension_sizes", std::vector<int>());
  _subregions = config_->get_or_default("subregions", std::vector<int>());
  init();
}
DataSource::~DataSource() {}

void DataSource::init() {
  std::cout << "DataSource initializing" << std::endl;
  // Not really sure what to do here...
}

void DataSource::populate(Mesh * const mesh_){
  if (mesh_->n_dimensions() != 2) {
    throw std::logic_error("Incorrect dimensionality for mesh type");
  }
  double *u0 = mesh_->get_u0();
  double *u1 = mesh_->get_u1();
  int cell_count = mesh_->get_cell_count();
  //N.B. zero is the 'former' matrix so this is where we place initial conds.
  size_t subregion_count = _subregions.size() % 4;
  if (subregion_count > 0) {
    for (size_t i = 0; i < subregions; ++i) {
      double  xmin = *it++;
      double ymin = *it++;
      double xmax = *it++;
      double xmax = *it++;
      for (size_t y = 0; y < /*todo get mesh padded size vector. */) {
        for (size_t x = 0; x < /*ditto */ )  {
          if mesh->displacement_from_origin(0, y) > // method which, given a dimension and a cell number, computes offset from dimension origin

          if 

        }
      }
    }
  }




  for (

  std::fill(&u0[0], &u0[cell_count], 0);








  // Zero initialize u1 - not strictly necessary
  std::fill(&u1[0], &u1[cell_count], 0);
};


