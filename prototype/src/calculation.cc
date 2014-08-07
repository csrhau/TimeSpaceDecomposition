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
  const double dx = _mesh->get_del_x();
  const double dy = _mesh->get_del_y();
  const double rx = dt_ / (dx * dx);
  const double ry = dt_ / (dy * dy);
  const int x_span = _mesh->get_node_augmented_col_count();
  const int core_rows = _mesh->get_node_core_row_count();
  const int core_cols = _mesh->get_node_core_col_count();
  const int i_offset = _mesh->get_current_row_offset();
  const int j_offset = _mesh->get_current_col_offset();
  for (int i = i_offset; i < core_rows + i_offset; ++i) {
    for (int j = j_offset; j < core_cols + j_offset; ++j) {
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
