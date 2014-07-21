#include "mesh.h"
#include "config_file.h"

#include <vector>
#include <iostream>
#include <algorithm>

Mesh::Mesh(const ConfigFile * const config_) : _config(config_) {
  _debug = config_->get_or_default("debug", false);
  _n_dimensions = config_->get_or_default("n_dimensions", 2);
  _dim_sizes = config_->get_or_default("dimension_sizes",
                                     std::vector<int>()); // Default will cause crash. Should be validated in main

  _padded_sizes = std::vector<int>();
  std::vector<int>::const_iterator it = _dim_sizes.begin();
  std::vector<int>::const_iterator itEnd = _dim_sizes.end();
  for (; it != itEnd; ++it) {
    _padded_sizes.push_back(*it + 2);
  }
  init();
}
Mesh::~Mesh() {}

void Mesh::init() {
  std::cout << "Mesh initializing" << std::endl;
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
    std::cout << "Constructed mesh of " << _cell_count << " cells" << std::endl;
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



