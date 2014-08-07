#include "static_mesh.h"

#include <mpi.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "tools-inl.h"
#include "config_file.h"

StaticMesh::StaticMesh(const ConfigFile& config_, 
           MPI_Comm cart_comm_,
           const std::vector<int>& dim_nodes_) : DistributedMesh(config_, 
                                                                cart_comm_,
                                                                dim_nodes_) {
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

StaticMesh::~StaticMesh() {
  delete[] _u0;
  delete[] _u1;
}

void StaticMesh::reflect_boundary(int boundary_) {
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

void StaticMesh::advance() {
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

void StaticMesh::exchange_boundaries() {
  // Use a very simple | 0 | 1 | 0 | 1 | scheme
  // 0 sends right, 1 sends left, then flip
  const int x_span = get_node_augmented_col_count();
  const int horizontal_cells = get_node_core_col_count();
  int paircount = 0;
  MPI_Request send_request[4]; // Has to be present but are not consulted
  MPI_Request recv_request[4];
  MPI_Status  recv_status[4];
  // TOP 
  if (has_top_neighbour()){
    const int i = 0;
    const int j = 1;
    MPI_Isend(&_u1[(i + 1) * x_span + j], horizontal_cells, MPI_DOUBLE, get_neighbour_rank(TOP), TOP, _cart_comm, &send_request[paircount]);
    MPI_Irecv(&_u1[i * x_span + j], horizontal_cells, MPI_DOUBLE, get_neighbour_rank(TOP), BOTTOM, _cart_comm, &recv_request[paircount]);
    ++paircount;
  }
  // LEFT
  if (has_left_neighbour()) {
    const int i = 1;
    const int j = 0;
    // irecv to i, j; isend from i, j+1
    MPI_Isend(&_u1[i * x_span + (j + 1)], 1, _col_type, get_neighbour_rank(LEFT), LEFT, _cart_comm, &send_request[paircount]);
    MPI_Irecv(&_u1[i * x_span + j], 1, _col_type, get_neighbour_rank(LEFT), RIGHT, _cart_comm, &recv_request[paircount]);
    ++paircount;
  }
  // BOTTOM
  if (has_bottom_neighbour()){
    const int i = get_node_core_row_count();
    const int j = 1;
    MPI_Isend(&_u1[i * x_span + j], horizontal_cells, MPI_DOUBLE, get_neighbour_rank(BOTTOM), BOTTOM, _cart_comm, &send_request[paircount]);
    MPI_Irecv(&_u1[(i + 1) * x_span + j], horizontal_cells, MPI_DOUBLE, get_neighbour_rank(BOTTOM), TOP, _cart_comm, &recv_request[paircount]);
    ++paircount;
  }
  // RIGHT
  if (has_right_neighbour()) {
    const int i = 1; 
    const int j = get_node_core_col_count();
    // irecv to i, j+1; isend from i, j
    MPI_Isend(&_u1[i * x_span + j], 1, _col_type, get_neighbour_rank(RIGHT), RIGHT, _cart_comm, &send_request[paircount]);
    MPI_Irecv(&_u1[i * x_span + (j + 1)], 1, _col_type, get_neighbour_rank(RIGHT), LEFT, _cart_comm, &recv_request[paircount]);
    ++paircount;
  }

  for (int req = 0; req < paircount; ++req) {
    MPI_Request_free(&send_request[req]);
  }
  MPI_Waitall(paircount, recv_request, recv_status);
}

double * StaticMesh::get_u0() { return _u0; }
double * StaticMesh::get_u1() { return _u1; }
int StaticMesh::get_node_core_row_count() const { return _node_core_row_count; }
int StaticMesh::get_node_core_col_count() const { return _node_core_col_count; }
int StaticMesh::get_node_augmented_row_count() const { return _node_augmented_row_count; }
int StaticMesh::get_node_augmented_col_count() const { return _node_augmented_col_count; }
int StaticMesh::get_node_core_cell_count() const {
  return get_node_core_row_count() * get_node_core_col_count();
}

int StaticMesh::get_node_augmented_cell_count() const {
  return get_node_augmented_row_count() * get_node_augmented_col_count();
}

int StaticMesh::get_current_row_offset() const {
  return 1;
}

int StaticMesh::get_current_col_offset() const {
  return 1;
}

int StaticMesh::get_previous_row_offset() const {
  return 1;
}

int StaticMesh::get_previous_col_offset() const {
  return 1;
}

// NOTE!! be very careful with augmented_ vs core in the following functions as
// core are logically 1-padded around the boundary
// core assume that 0 is the logical first column..
double StaticMesh::get_core_row_y(int row_) const {
  return _core_origin_y + row_ * get_del_y();
}
double StaticMesh::get_core_col_x(int col_) const {
  return _core_origin_x + col_ * get_del_x();
}
// NOTE - could be re-written to use augmented_origin and not subtract 1 from row
double StaticMesh::get_y_coord(int row_) const {
  return _core_origin_y + (row_ - 1) * get_del_y();
}
double StaticMesh::get_x_coord(int col_) const {
  return _core_origin_x + (col_ - 1) * get_del_x();
}
