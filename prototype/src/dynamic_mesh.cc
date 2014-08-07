#include "dynamic_mesh.h"

#include <mpi.h>
#include <iostream>

#include "tools-inl.h"
#include "config_file.h"

DynamicMesh::DynamicMesh(const ConfigFile& config_,
                         MPI_Comm cart_comm_,
                         const std::vector<int>& dim_nodes_) : DistributedMesh(config_,
                                                                               cart_comm_,
                                                                               dim_nodes_),
                                                               _prograde(true) {
  const int memory_elements = get_node_augmented_cell_count();
  _u0 = new double[memory_elements];
  _u1 = new double[memory_elements];
  // Create MPI Datatypes
}

DynamicMesh::~DynamicMesh() {
  delete[] _u0;
  delete[] _u1;
}


void DynamicMesh::exchange_boundaries() {
  // Plan: Send padded rows, 2 down on prograde, 2 up on retrograde
  MPI_Request send_request[4]; // Has to be present but are not consulted
  MPI_Request recv_request[4];
  MPI_Status  recv_status[4];
  int send_count = 0;
  int recv_count = 0;
  // This will have to change when we're dealing in 2d properly
  const int x_span = get_node_augmented_col_count();
  if (_prograde) {
    // SEND
    if (has_bottom_neighbour()) {
      const int i = get_node_augmented_row_count() - 3; // one for last ghost
                                                        // and 2 to send 2
      const int j = 0;
      MPI_Isend(&_u1[i * x_span + j],
                2 * x_span,
                MPI_DOUBLE,
                get_neighbour_rank(BOTTOM),
                BOTTOM,
                _cart_comm,
                &send_request[send_count++]);
    }
    // RECEIVE
    if (has_top_neighbour()) {
      const int i = 0;
      const int j = 0;
      MPI_Irecv(&_u1[i * x_span + j],
                2 * x_span,
                MPI_DOUBLE,
                get_neighbour_rank(TOP),
                BOTTOM,
                _cart_comm,
                &recv_request[recv_count++]);
    }
  } else { /* retrograde */
    // SEND
    if (has_top_neighbour()) {
      const int i = 1; // Start of core rows...
      const int j = 0;
      MPI_Isend(&_u1[i * x_span + j],
                2 * x_span,
                MPI_DOUBLE,
                get_neighbour_rank(TOP),
                TOP,
                _cart_comm,
                &send_request[send_count++]);
    }
    // RECEIVE
    if (has_bottom_neighbour()) {
      const int i = get_node_augmented_row_count() - 2;
      const int j = 0;
      MPI_Irecv(&_u1[i * x_span + j],
                2 * x_span,
                MPI_DOUBLE,
                get_neighbour_rank(BOTTOM),
                TOP,
                _cart_comm,
                &recv_request[recv_count++]);
    }
  }
  for (int req = 0; req < send_count; ++req) {
    MPI_Request_free(&send_request[req]);
  }
  MPI_Waitall(recv_count, recv_request, recv_status);
}

void DynamicMesh::advance() {
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
  // Toggle step
  _prograde = !_prograde;
  // Now we've finished updating u1, we can swap it to u0
  std::swap(_u0, _u1);
}

void DynamicMesh::reflect_boundary(int boundary_) {
  const int x_span = get_node_augmented_col_count();
  const int core_rows = get_node_core_row_count();
  const int core_cols = get_node_core_col_count();
  const int row_offset = get_current_row_offset();
  const int col_offset = get_current_col_offset();
  switch (boundary_) {
    case (TOP): {
      int i = row_offset;
      for (int j = col_offset;
          j < core_cols + col_offset;
          ++j) {
        const int top = (i - 1) * x_span + j;
        const int center = i * x_span + j;
        _u1[top] = _u1[center];
      }
    }; break;
    case (BOTTOM): {
      int i = core_rows + row_offset - 1;
      for (int j = col_offset;
          j < core_cols + col_offset;
          ++j) {
        const int bottom = (i + 1) * x_span + j;
        const int center = i * x_span + j;
        _u1[bottom] = _u1[center];
      }
    }; break;
    case (LEFT): {
      int j = col_offset;
      for (int i = row_offset; i < core_rows + row_offset; ++i) {
        const int left = i * x_span + (j - 1);
        const int center = i * x_span + j;
        _u1[left] = _u1[center];
      }
    }; break;
    case (RIGHT): {
      int j = core_cols + col_offset - 1;
      for (int i = row_offset; i < core_rows + row_offset; ++i) {
        const int right = i * x_span + (j + 1);
        const int center = i * x_span + j;
        _u1[right] = _u1[center];
      }
    }; break;
  }
}

double * DynamicMesh::get_u0() { return _u0; }
double * DynamicMesh::get_u1() { return _u1; }

int DynamicMesh::get_current_row_offset() const {
  if (_prograde && has_top_neighbour()) {
    return 2;
  } else {
    return 1;
  }
}

int DynamicMesh::get_current_col_offset() const {
  if (_prograde && has_left_neighbour()) {
    return 2;
  } else {
    return 1;
  }
}

// Previous versions simply negate prograde
int DynamicMesh::get_previous_row_offset() const {
  if (!_prograde && has_top_neighbour()) {
    return 2;
  } else {
    return 1;
  }
}

int DynamicMesh::get_previous_col_offset() const {
  if (!_prograde && has_left_neighbour()) {
    return 2;
  } else {
    return 1;
  }
}

// Varies w.r.t. prograde/retrograde
double DynamicMesh::get_core_origin_x() const {
  return get_del_x() * calculate_local_offset(get_node_col(),
                                              get_horizontal_nodes_count(),
                                              get_world_core_col_count(),
                                              _prograde);
}
double DynamicMesh::get_core_origin_y() const {
  return get_del_y() * calculate_local_offset(get_node_row(),
                                              get_vertical_nodes_count(),
                                              get_world_core_row_count(),
                                              _prograde);
}

// Invariant w.r.t. prograde/retrograde - deals with 'outer/padded' notation, i.e. just raw indexes
double DynamicMesh::get_x_coord(int col_) const {
  return get_core_origin_x() + (col_ - 1) * get_del_x();
}

// Invariant w.r.t. prograde/retrograde - deals with 'outer/padded' notation, i.e. just raw indexes
double DynamicMesh::get_y_coord(int row_) const {
  return get_core_origin_y() + (row_ - 1) * get_del_y();
}

// Varies given prograde/retrograde
int DynamicMesh::get_node_core_row_count() const {
  return calculate_local_span(get_node_row(),
                              get_vertical_nodes_count(),
                              get_world_core_row_count(),
                              _prograde);
}
// Varies given prograde/retrograde
int DynamicMesh::get_node_core_col_count() const {
  return calculate_local_span(get_node_col(),
                              get_horizontal_nodes_count(),
                              get_world_core_col_count(),
                              _prograde);
}
// Varies given prograde/retrograde
int DynamicMesh::get_node_core_cell_count() const {
  return get_node_core_row_count()
       * get_node_core_col_count();
}
// constant, invariant w.r.t. prograde/retrograde
int DynamicMesh::get_node_augmented_row_count() const {
  const int border_rows = 2;
  const int unaugmented_rows = calculate_local_span(get_node_row(),
                                                    get_vertical_nodes_count(),
                                                    get_world_core_row_count());
  const int additional_rows = has_top_neighbour() ? 1 : 0;
  return unaugmented_rows + border_rows + additional_rows;
}
// constant, invariant w.r.t. prograde/retrograde
int DynamicMesh::get_node_augmented_col_count() const {
  const int border_cols = 2;
  const int unaugmented_cols = calculate_local_span(get_node_col(),
                                                    get_horizontal_nodes_count(),
                                                    get_world_core_col_count());
  const int additional_cols = has_left_neighbour() ? 1 : 0;
  return unaugmented_cols + border_cols + additional_cols;
}
// constant, invariant w.r.t. prograde/retrograde
int DynamicMesh::get_node_augmented_cell_count() const {
  return get_node_augmented_row_count()
       * get_node_augmented_col_count();
}
