#include "mesh.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "config_file.h"
#include "tools.h"

Mesh::Mesh(const ConfigFile * const config_) : _config(config_) {
  parse_config();
  _padded_sizes = std::vector<int>();
  std::vector<int>::const_iterator it = _dim_cells.begin();
  std::vector<int>::const_iterator itEnd = _dim_cells.end();
  for (; it != itEnd; ++it) {
    _from_indexes.push_back(1);
    _to_indexes.push_back(*it + 1);
    _padded_sizes.push_back(*it + 2);
  }
  init();
}

Mesh::~Mesh() {
  delete[]  _u0;
  delete[] _u1;
}

void Mesh::parse_config() {
  _debug = _config->get_or_default("debug", false);
  _n_dimensions = _config->get_or_default("n_dimensions", 0);
  _min_coords = _config->get_or_default("min_coordinates", std::vector<double>());
  _max_coords = _config->get_or_default("max_coordinates", std::vector<double>());
  _dim_cells = _config->get_or_default("dimension_cells",
                                     std::vector<int>()); // Default will cause crash. Should be validated in main
}

void Mesh::init() {
  for (int dim = 0; dim < _n_dimensions; ++dim) {
    std::vector<double> offsets;
    double delta = (_max_coords[dim] - _min_coords[dim]) / _dim_cells[dim];
    _dimension_deltas.push_back(delta);
    for (int cell = 0; cell < _padded_sizes[dim]; ++cell) {
      double offset = _min_coords[dim] + delta * (cell - 1);
      // TODO!!!!! the -1 takes into account lhs padding. Not guaranteed
      // to always be there
      offsets.push_back (offset);
    }
    _cell_offsets.push_back(offsets);
  }
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

void Mesh::update_boundaries(){
  for (int boundary = 0; boundary < 4; ++boundary) {
    reflect_boundary(boundary);
  }
} // Should be virtual

void Mesh::reflect_boundary(int boundary_) {
  // N.B. use u1 as we're in the current timestep
  double *u1 = get_u1();
  int x_min = get_from_index(DIM_X);
  int y_min = get_from_index(DIM_Y);
  int x_max = get_to_index(DIM_X);
  int y_max = get_to_index(DIM_Y);
  int x_span = get_dimension_span(DIM_X);
  switch(boundary_) {
    case 0: { /* bottom */
      for (int i = x_min; i < x_max; ++i) {
          int bottom = POLY2(i, y_max, x_span);
          int center = POLY2(i, y_max-1, x_span);
          u1[bottom] = u1[center];
      }
    } break;
    case 1: { /* left */
      for (int j = y_min; j < y_max; ++j) {
        int left = POLY2(x_min - 1, j, x_span);
        int center = POLY2(x_min, j, x_span);
        u1[left] = u1[center];
      }
    } break;
    case 2: { /* top */
      for (int i = x_min; i < x_max; ++i) {
          int top = POLY2(i, y_min - 1, x_span);
          int center = POLY2(i, y_min, x_span);
          u1[top] = u1[center];
      }
    } break;
    case 3: { /* right */
      for (int j = y_min; j < y_max; ++j) {
        int right = POLY2(x_max, j, x_span);
        int center = POLY2(x_max - 1, j, x_span);
        u1[right] = u1[center];
      }
    } break;
    default: {
      throw std::domain_error("Error in reflectBoundaries(): unknown boundary "
                              + boundary_);
    }
  }
}

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

int Mesh::get_from_index(int dim_) const {
  return _from_indexes.at(dim_);
}

int Mesh::get_to_index(int dim_) const {
  return _to_indexes.at(dim_);
}

double Mesh::get_dimension_delta(int dim_) const {
  return  _dimension_deltas.at(dim_);
}

int Mesh::get_internal_cells(int dim_) const {
  return _dim_cells.at(dim_);
}

// Number of cells for a "row" of the space - n.b. includes padding
int Mesh::get_dimension_span(int dim_) const {
  return _padded_sizes.at(dim_);
}
