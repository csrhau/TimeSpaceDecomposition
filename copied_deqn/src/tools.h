#ifndef TOOLS_H
#define TOOLS_H

#define POLY2(i, j, ispan) ((i) + ((j)  * (ispan)))

#define DIM_X 0
#define DIM_Y 1
#define DIM_Z 2

namespace Tools {

static inline double mesh_sum_2d(Mesh* mesh_) {
  double *u0 = mesh_->get_u0();
  int x_min = mesh_->get_from_index(DIM_X);
  int y_min = mesh_->get_from_index(DIM_Y);
  int x_max = mesh_->get_to_index(DIM_X);
  int y_max = mesh_->get_to_index(DIM_Y);
  int x_span = mesh_->get_dimension_span(DIM_X);
  double sum = 0;
  for (int j = y_min; j < y_max; ++j) {
    for (int i = x_min; i < x_max; ++i) {
      int center = POLY2(i, j, x_span);
      sum += u0[center];
    }
  }
  return sum;
}
} // namespace Tools
#endif
