#include "calculation.h"

#include "config_file.h"
#include "mesh.h"

Calculation::Calculation(const ConfigFile& config_, Mesh *mesh_)
                        : _config(config_), _mesh(mesh_) {}

Calculation::~Calculation() {}

void Calculation::step(double dt_) {
  diffuse(dt_);
}

void Calculation::diffuse(double dt_) {
  double *u0 = _mesh->get_u0();
  double *u1 = _mesh->get_u1();
  double dx = _mesh->get_del_x();
  double dy = _mesh->get_del_y();
  double rx = dt_ / (dx * dx);
  double ry = dt_ / (dy * dy);
  int x_span = _mesh->get_node_outer_cols();
  for (int i = 1; i < _mesh->get_node_inner_rows() + 1; ++i) {
    for (int j = 1; j < _mesh->get_node_inner_cols() + 1; ++j) {
      const int center = i * x_span + j;
      const int top = (i - 1) * x_span + j;
      const int bottom = (i + 1) * x_span + j;
      const int left = i * x_span + (j - 1);
      const int right = i * x_span + (j + 1);
      u1[center] = (1.0 - 2.0*rx - 2.0*ry) *u0[center] + rx * u0[left]
                         + rx * u0[right] + ry * u0[top] + ry * u0[bottom];
    }
  }
}
