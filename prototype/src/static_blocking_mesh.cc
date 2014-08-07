#include "static_blocking_mesh.h"

#include <mpi.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "tools-inl.h"
#include "config_file.h"

StaticBlockingMesh::StaticBlockingMesh(const ConfigFile& config_, 
           MPI_Comm cart_comm_,
           const std::vector<int>& dim_nodes_) : DistributedMesh(config_,
                                                                 cart_comm_,
                                                                 dim_nodes_) {
  // Calculate our simulation domain
  // Space not including ghost cells and boundary padding

  _node_core_row_count = calculate_local_span(get_node_row(),
                                          get_vertical_nodes_count(),
                                          get_world_core_row_count());
  _node_core_col_count = calculate_local_span(get_node_col(),
                                          get_horizontal_nodes_count(),
                                          get_world_core_col_count());
  // TODO: If we ever want non-zero origins, we need to add the offset here
  _core_origin_y = get_del_y() * calculate_local_offset(get_node_row(),
                                                         get_vertical_nodes_count(),
                                                         get_world_core_row_count());

  _core_origin_x = get_del_x() * calculate_local_offset(get_node_col(),
                                                        get_horizontal_nodes_count(),
                                                        get_world_core_col_count());


  _node_augmented_row_count = _node_core_row_count + 2;
  _node_augmented_col_count = _node_core_col_count + 2;
  _u0 = new double[_node_augmented_row_count * _node_augmented_col_count];
  _u1 = new double[_node_augmented_row_count * _node_augmented_col_count];
  // Create MPI vector datatype to deal with column sending.
  MPI_Type_vector(_node_core_row_count, // # column height
                  1,                // 1 column only
                  _node_augmented_col_count, // x dimension span
                  MPI_DOUBLE,
                  &_col_type);
  MPI_Type_commit(&_col_type);
}

StaticBlockingMesh::~StaticBlockingMesh() {
  delete[] _u0;
  delete[] _u1;
}

void StaticBlockingMesh::reflect_boundary(int boundary_) {
  // n.b. use u1 as we're in the current timestep
  int x_span = get_node_augmented_col_count();
  switch (boundary_) {
    case (TOP): {
      int i = 1;
      for (int j = 1; j < get_node_core_col_count() + 1; ++j) {
        const int top = (i - 1) * x_span + j;
        const int center = i * x_span + j;
        _u1[top] = _u1[center];
      } 
    } break;
    case (BOTTOM): {
      int i = get_node_core_row_count(); // - n.b. this includes bound
      for (int j = 1; j < get_node_core_col_count() + 1; ++j) {
        const int bottom = (i + 1) * x_span + j;
        const int center = i * x_span + j;
        _u1[bottom] = _u1[center];
      } 
    } break;
    case (LEFT): {
      int j = 1; 
      for (int i = 1; i < get_node_core_row_count() + 1; ++i) {
        const int left = i * x_span + (j - 1);
        const int center = i * x_span + j;
        _u1[left] = _u1[center];
      }
    } break;
    case (RIGHT): {
      int j = get_node_core_col_count(); // again includes bound
      for (int i = 1; i < get_node_core_row_count() + 1; ++i) {
        const int right = i * x_span + (j + 1);
        const int center = i * x_span + j;
        _u1[right] = _u1[center];
      }
    } break;
  }
}

void StaticBlockingMesh::advance() {
  if (!has_top_neighbour()) {
    reflect_boundary(TOP);
  }
  if (!has_bottom_neighbour()) {
    reflect_boundary(BOTTOM);
  }
  if (!has_left_neighbour()) {
    reflect_boundary(LEFT);
  }
  if (!has_right_neighbour()) {
    reflect_boundary(RIGHT);
  }
  exchange_boundaries();
  // Now we've finished updating u1, we can swap it to u0
  std::swap(_u0, _u1);
}

void StaticBlockingMesh::exchange_boundaries() {
  // Plan: Use sendrecv as follows:
  // Colwise, odd up, even down          - TAGGED TOP
  // Colwise, odd down, even up          - TAGGED BOTTOM
  // Rowwise, odd left, even right       - TAGGED LEFT
  // Rowwise, odd right, even left       - TAGGED RIGHT
  MPI_Status status;
  const int x_span = get_node_augmented_col_count();
  const int horizontal_cells = get_node_core_col_count();
  const int vertical_cells = get_node_core_row_count();
  const bool odd_row = get_node_row() & 1;
  const bool odd_col = get_node_col() & 1;
  double * const up_sendbuf    = &_u1[1 * x_span + 1];
  double * const down_sendbuf  = &_u1[vertical_cells * x_span + 1];
  double * const left_sendbuf  = up_sendbuf;
  double * const right_sendbuf = &_u1[x_span + horizontal_cells];
  double * const up_recvbuf    = &_u1[1];
  double * const down_recvbuf  = &_u1[(vertical_cells + 1) * x_span + 1];
  double * const left_recvbuf  = &_u1[x_span];
  double * const right_recvbuf = &_u1[x_span + horizontal_cells + 1]; 

  // STEP 1: Odd rows up, evens down. Note - odds always have top neighbour
  if ((odd_row && has_top_neighbour()) || (!odd_row && has_bottom_neighbour())) {
    MPI_Sendrecv(odd_row ? up_sendbuf : down_sendbuf,
                 horizontal_cells,
                 MPI_DOUBLE,
                 odd_row ? get_neighbour_rank(TOP) : get_neighbour_rank(BOTTOM),
                 TOP,
                 odd_row ? up_recvbuf : down_recvbuf,
                 horizontal_cells,
                 MPI_DOUBLE,
                 odd_row ? get_neighbour_rank(TOP) : get_neighbour_rank(BOTTOM),
                 TOP,
                 _cart_comm,
                 &status);

  }
  // STEP 2: Odd rows down, evens up
  if ((odd_row && has_bottom_neighbour()) || (!odd_row && has_top_neighbour())) {
    MPI_Sendrecv(odd_row ? down_sendbuf : up_sendbuf,
                 horizontal_cells,
                 MPI_DOUBLE,
                 odd_row ? get_neighbour_rank(BOTTOM) : get_neighbour_rank(TOP),
                 BOTTOM,
                 odd_row ? down_recvbuf : up_recvbuf,
                 horizontal_cells,
                 MPI_DOUBLE,
                 odd_row ? get_neighbour_rank(BOTTOM) : get_neighbour_rank(TOP),
                 BOTTOM,
                 _cart_comm,
                 &status);
  }
  // STEP 3: Odd cols left, evens right
  if ((odd_col && has_left_neighbour()) || (!odd_col && has_right_neighbour())) {
  MPI_Sendrecv(odd_col ? left_sendbuf : right_sendbuf,
               1,
               _col_type,
               odd_col ? get_neighbour_rank(LEFT) : get_neighbour_rank(RIGHT),
               LEFT,
               odd_col ? left_recvbuf : right_recvbuf,
               1,
               _col_type,
               odd_col ? get_neighbour_rank(LEFT) : get_neighbour_rank(RIGHT),
               LEFT,
               _cart_comm,
               &status);
  }
  // STEP 3: Odd cols right, evens left
  if ((odd_col && has_right_neighbour()) || (!odd_col && has_left_neighbour())) {
  MPI_Sendrecv(odd_col ? right_sendbuf : left_sendbuf,
               1,
               _col_type,
               odd_col ? get_neighbour_rank(RIGHT): get_neighbour_rank(LEFT),
               RIGHT,
               odd_col ? right_recvbuf : left_recvbuf,
               1,
               _col_type,
               odd_col ? get_neighbour_rank(RIGHT): get_neighbour_rank(LEFT),
               RIGHT,
               _cart_comm,
               &status);
  }
}

double * StaticBlockingMesh::get_u0() { return _u0; }
double * StaticBlockingMesh::get_u1() { return _u1; }
int StaticBlockingMesh::get_node_core_row_count() const { return _node_core_row_count; }
int StaticBlockingMesh::get_node_core_col_count() const { return _node_core_col_count; }
int StaticBlockingMesh::get_node_augmented_row_count() const { return _node_augmented_row_count; }
int StaticBlockingMesh::get_node_augmented_col_count() const { return _node_augmented_col_count; }
int StaticBlockingMesh::get_node_core_cell_count() const {
  return get_node_core_row_count() * get_node_core_col_count();
}
int StaticBlockingMesh::get_node_augmented_cell_count() const {
  return get_node_augmented_row_count() * get_node_augmented_col_count();
}

int StaticBlockingMesh::get_current_row_offset() const {
  return 1;
}

int StaticBlockingMesh::get_current_col_offset() const {
  return 1;
}

int StaticBlockingMesh::get_previous_row_offset() const {
  return 1;
}

int StaticBlockingMesh::get_previous_col_offset() const {
  return 1;
}


// NOTE!! be very careful with augmented_ vs core in the following functions as
// core are logically 1-padded around the boundary
// core assume that 0 is the logical first column..
double StaticBlockingMesh::get_core_row_y(int row_) const {
  return _core_origin_y + row_ * get_del_y();
}
double StaticBlockingMesh::get_core_col_x(int col_) const {
  return _core_origin_x + col_ * get_del_x();
}
// NOTE - could be re-written to use augmented_origin and not subtract 1 from row
double StaticBlockingMesh::get_y_coord(int row_) const {
  return _core_origin_y + (row_ - 1) * get_del_y();
}
double StaticBlockingMesh::get_x_coord(int col_) const {
  return _core_origin_x + (col_ - 1) * get_del_x();
}


