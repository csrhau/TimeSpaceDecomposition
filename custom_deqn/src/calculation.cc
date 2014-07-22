#include "calculation.h"

#include "config_file.h"
#include "mesh.h"

#include <iostream>

#define POLY2(i, j, ispan) ((i) + ((j)  * (ispan)))

Calculation::Calculation(const ConfigFile * const config_,
                         Mesh * const mesh_) : _config(config_),
                                               _mesh(mesh_) {
  _n_dimensions = _config->get_or_default("n_dimensions", 2);
}

Calculation::~Calculation() {}

void Calculation::init() {}

void Calculation::step(double dt_) {
  diffuse(dt_);
}

// Should be virtual
void Calculation::diffuse(double dt_) {
  double *u0 = _mesh->get_u0();
  double *u1 = _mesh->get_u1();
  int x_min = _mesh->get_from_index(0);
  int y_min = _mesh->get_from_index(1);
  int x_max = _mesh->get_to_index(0);
  int y_max = _mesh->get_to_index(1);
  int x_span = _mesh->get_dimension_span(0);
  double dx = _mesh->get_dimension_delta(0);
  double dy = _mesh->get_dimension_delta(1);
  double rx = dt_ / (dx * dx);
  double ry = dt_ / (dy * dy);
  for (int j = y_min; j < y_max; ++j) {
    for (int i = x_min; i < x_max; ++i) {
      int top    = POLY2(i, j-1, x_span);
      int left   = POLY2(i-1, j, x_span);
      int center = POLY2(i, j, x_span);
      int right  = POLY2(i+1, j, x_span);
      int bottom = POLY2(i, j+1, x_span);
      u1[center] = (1.0 - 2.0*rx - 2.0*ry) *u0[center] + rx * u0[left] 
                 + rx * u0[right] + ry * u0[top] + ry * u0[bottom];
    }
  }
}

int  Calculation::n_dimensions() const { 
  return _n_dimensions;
}
