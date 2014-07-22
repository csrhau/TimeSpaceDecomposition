#include "mesh.h"
#include "config_file.h"

#include <vector>
#include <iostream>
#include <algorithm>

Mesh::Mesh(const ConfigFile * const config_) : _config(config_) {
  _debug = config_->get_or_default("debug", false);
  _n_dimensions = config_->get_or_default("n_dimensions", 0);
  _min_coords = config_->get_or_default("min_coordinates", std::vector<double>());
  _max_coords = config_->get_or_default("max_coordinates", std::vector<double>());
  _dim_cells = config_->get_or_default("dimension_cells",
                                     std::vector<int>()); // Default will cause crash. Should be validated in main
  _padded_sizes = std::vector<int>();
  std::vector<int>::const_iterator it = _dim_cells.begin();
  std::vector<int>::const_iterator itEnd = _dim_cells.end();
  for (; it != itEnd; ++it) {
    _padded_sizes.push_back(*it + 2);
  }
  init();
}
Mesh::~Mesh() {}

void Mesh::init() {
  for (int dim = 0; dim < _n_dimensions; ++dim) {
    std::vector<double> offsets;
    double delta = (_max_coords[dim] - _min_coords[dim]) / _dim_cells[dim];
    for (int cell = 0; cell < _padded_sizes[dim]; ++cell) {
      double offset = _min_coords[dim] + delta * (cell - 1);
      // TODO!!!!! the -1 takes into account lhs padding. Not guaranteed
      // to always be there
      offsets.push_back (offset);
    }
    _cell_offsets.push_back(offsets);
  } 
  // TODO do cell_offsets [0 = x, 1 = y and so on]
  // Calculate total memory requirements for each space
  std::vector<int>::const_iterator it = _padded_sizes.begin();
  std::vector<int>::const_iterator itEnd = _padded_sizes.end();
  _cell_count = *it;
  for (++it; it != itEnd; ++it) {
    _cell_count *= (*it);
  }
  _u0 = new double[_cell_count];
  _u1 = new double[_cell_count];
  if (_debug) {
    std::cout << "Allocated mesh of " << _cell_count << " cells" << std::endl;
  }
}

void Mesh::step() {
  std::swap(_u0, _u1);
} // Should be virtual
int Mesh::n_dimensions() const {
  return _n_dimensions;
}

double * Mesh::get_u0() {
  return _u0;
}

double * Mesh::get_u1() {
  return _u1;
}

int Mesh::get_cell_count() {
  return _cell_count;
}

double Mesh::displacement(int dimension, int cell) const {
  // Throws std::out_of_range exception.
  return _cell_offsets.at(dimension).at(cell);
}

const std::vector<int>& Mesh::get_padded_sizes() const {
  return _padded_sizes;
}
