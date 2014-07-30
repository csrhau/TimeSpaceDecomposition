#include "mesh.h"

#include <mpi.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "tools-inl.h"
#include "config_file.h"

namespace {
  enum Boundary {
    TOP = 1, 
    BOTTOM = 2, 
    LEFT = 3,
    RIGHT = 4,
  };
}

Mesh::Mesh(const ConfigFile& config_, 
           MPI_Comm cart_comm_,
           const std::vector<int>& dim_nodes_) : _config(config_), 
                                                 _cart_comm(cart_comm_),
                                                 _dim_nodes(dim_nodes_) {
  MPI_Comm_rank(_cart_comm, &_cart_rank);
  // Calculate our simulation domain
  // Space not including ghost cells and boundary padding
  std::vector<int> inner_dimensions = _config.get_or_default("logical_dimensions",
                                                             std::vector<int>());
  _world_inner_rows = inner_dimensions.at(0);
  _world_inner_cols = inner_dimensions.at(1);
  std::vector<double> physical_dimensions = _config.get_or_default("physical_dimensions",
                                                             std::vector<double>());
  _world_height = physical_dimensions.at(0);
  _world_width = physical_dimensions.at(1);

  _cart_coords.resize(2, 0); 
  MPI_Cart_coords(_cart_comm, _cart_rank, 2, &_cart_coords[0]);

  _node_inner_rows = calculate_local_span(_cart_coords.at(0),
                                          _dim_nodes.at(0),
                                          inner_dimensions.at(0));
  _node_inner_cols = calculate_local_span(_cart_coords.at(1),
                                          _dim_nodes.at(1),
                                          inner_dimensions.at(1));
  // TODO: If we ever want non-zero origins, we need to add the offset here
  _inner_origin_y = get_del_y() * calculate_local_offset(_cart_coords.at(0),
                                                         _dim_nodes.at(0),
                                                         inner_dimensions.at(0));

  _inner_origin_x = get_del_x() * calculate_local_offset(_cart_coords.at(1),
                                                         _dim_nodes.at(1),
                                                         inner_dimensions.at(1));


  _node_outer_rows = _node_inner_rows + 2;
  _node_outer_cols = _node_inner_cols + 2;
  _u0 = new double[_node_outer_rows * _node_outer_cols];
  _u1 = new double[_node_outer_rows * _node_outer_cols];
  // Create MPI vector datatype to deal with column sending.
  MPI_Type_vector(_node_inner_rows, // # column height
                  1,                // 1 column only
                  _node_outer_cols, // x dimension span
                  MPI_DOUBLE,
                  &_col_type);
  MPI_Type_commit(&_col_type);

  // Compute ranks of neighbours
  int coords[2];
  if (has_top_neighbour()) {
    coords[0] = _cart_coords[0] - 1;
    coords[1] = _cart_coords[1];
    MPI_Cart_rank(_cart_comm, coords, &_top_rank_or_neg);
  }
  if (has_bottom_neighbour()) {
    coords[0] = _cart_coords[0] + 1;
    coords[1] = _cart_coords[1];
    MPI_Cart_rank(_cart_comm, coords, &_bottom_rank_or_neg);
  }
  if (has_left_neighbour()) {
    coords[0] = _cart_coords[0];
    coords[1] = _cart_coords[1] - 1;
    MPI_Cart_rank(_cart_comm, coords, &_left_rank_or_neg);
  }
  if (has_right_neighbour()) {
    coords[0] = _cart_coords[0];
    coords[1] = _cart_coords[1] + 1;
    MPI_Cart_rank(_cart_comm, coords, &_right_rank_or_neg);
  }
}

Mesh::~Mesh() {
  delete[] _u0;
  delete[] _u1;
}

void Mesh::reflect_boundary(int boundary_) {
  // n.b. use u1 as we're in the current timestep
  int x_span = get_node_outer_cols();
  switch (boundary_) {
    case (TOP): {
      int i = 1;
      for (int j = 1; j < get_node_inner_cols() + 1; ++j) {
        const int top = (i - 1) * x_span + j;
        const int center = i * x_span + j;
        _u1[top] = _u1[center];
      } 
    } break;
    case (BOTTOM): {
      int i = get_node_inner_rows(); // - n.b. this includes bound
      for (int j = 1; j < get_node_inner_cols() + 1; ++j) {
        const int bottom = (i + 1) * x_span + j;
        const int center = i * x_span + j;
        _u1[bottom] = _u1[center];
      } 
    } break;
    case (LEFT): {
      int j = 1; 
      for (int i = 1; i < get_node_inner_rows() + 1; ++i) {
        const int left = i * x_span + (j - 1);
        const int center = i * x_span + j;
        _u1[left] = _u1[center];
      }
    } break;
    case (RIGHT): {
      int j = get_node_inner_cols(); // again includes bound
      for (int i = 1; i < get_node_inner_rows() + 1; ++i) {
        const int right = i * x_span + (j + 1);
        const int center = i * x_span + j;
        _u1[right] = _u1[center];
      }
    } break;
  }
}

void Mesh::advance() {
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

void Mesh::exchange_boundaries() {
  // Use a very simple | 0 | 1 | 0 | 1 | scheme
  // 0 sends right, 1 sends left, then flip
  const int x_span = get_node_outer_cols();
  const int horizontal_cells = get_node_inner_cols();
  int paircount = 0;
  MPI_Request send_request[4]; // Has to be present but are not consulted
  MPI_Request recv_request[4];
  MPI_Status  recv_status[4];
  if (has_top_neighbour()){
    const int i = 0;
    const int j = 1;
    MPI_Isend(&_u1[(i + 1) * x_span + j], horizontal_cells, MPI_DOUBLE, _top_rank_or_neg, TOP, _cart_comm, &send_request[paircount]);
    MPI_Irecv(&_u1[i * x_span + j], horizontal_cells, MPI_DOUBLE, _top_rank_or_neg, BOTTOM, _cart_comm, &recv_request[paircount]);
    ++paircount;
  }
  if (has_bottom_neighbour()){
    const int i = get_node_inner_rows();
    const int j = 1;
    MPI_Isend(&_u1[i * x_span + j], horizontal_cells, MPI_DOUBLE, _bottom_rank_or_neg, BOTTOM, _cart_comm, &send_request[paircount]);
    MPI_Irecv(&_u1[(i + 1) * x_span + j], horizontal_cells, MPI_DOUBLE, _bottom_rank_or_neg, TOP, _cart_comm, &recv_request[paircount]);
    ++paircount;
  }
  if (has_left_neighbour()) {
    const int i = 1;
    const int j = 0;
    // irecv to i, j; isend from i, j+1
    MPI_Isend(&_u1[i * x_span + (j + 1)], 1, _col_type, _left_rank_or_neg, LEFT, _cart_comm, &send_request[paircount]);
    MPI_Irecv(&_u1[i * x_span + j], 1, _col_type, _left_rank_or_neg, RIGHT, _cart_comm, &recv_request[paircount]);
    ++paircount;
  }
  if (has_right_neighbour()) {
    const int i = 1; 
    const int j = get_node_inner_rows();
    // irecv to i, j+1; isend from i, j
    MPI_Isend(&_u1[i * x_span + j], 1, _col_type, _right_rank_or_neg, RIGHT, _cart_comm, &send_request[paircount]);
    MPI_Irecv(&_u1[i * x_span + (j + 1)], 1, _col_type, _right_rank_or_neg, LEFT, _cart_comm, &recv_request[paircount]);
    ++paircount;
  }

  for (int req = 0; req < paircount; ++req) {
    MPI_Request_free(&send_request[req]);
  }
  MPI_Waitall(paircount, recv_request, recv_status);
}

double * Mesh::get_u0() { return _u0; }
double * Mesh::get_u1() { return _u1; }
int Mesh::get_node_inner_rows() const { return _node_inner_rows; }
int Mesh::get_node_inner_cols() const { return _node_inner_cols; }
int Mesh::get_node_outer_rows() const { return _node_outer_rows; }
int Mesh::get_node_outer_cols() const { return _node_outer_cols; }
int Mesh::get_node_inner_cell_count() const {
  return get_node_inner_rows() * get_node_inner_cols();
}
int Mesh::get_node_outer_cell_count() const {
  return get_node_outer_rows() * get_node_outer_cols();
}
double Mesh::get_world_inner_rows() const { return _world_inner_rows; }
double Mesh::get_world_inner_cols() const { return _world_inner_cols; }
double Mesh::get_del_y() const { return _world_height / _world_inner_rows; }
double Mesh::get_del_x() const { return _world_width / _world_inner_cols; }
// NOTE!! be very careful with outer_ vs inner in the following functions as
// inner are logically 1-padded around the boundary
// inner assume that 0 is the logical first column..
double Mesh::get_inner_row_y(int row_) const {
  return _inner_origin_y + row_ * get_del_y();
}
double Mesh::get_inner_col_x(int col_) const {
  return _inner_origin_x + col_ * get_del_x();
}
// NOTE - could be re-written to use outer_origin and not subtract 1 from row
double Mesh::get_outer_row_y(int row_) const {
  return _inner_origin_y + (row_ - 1) * get_del_y();
}
double Mesh::get_outer_col_x(int col_) const {
  return _inner_origin_x + (col_ - 1) * get_del_x();
}

bool Mesh::has_top_neighbour() const {
  // We have a top neighbour if we are not on the first row
  return _cart_coords[0] != 0; 
}

bool Mesh::has_bottom_neighbour() const {
  // We have a bottom neighbour if we are not on the last row
  return _cart_coords[0] != _dim_nodes[0] - 1; 
}

bool Mesh::has_left_neighbour() const {
  // We have a left neighbour if we are not on the first column
  return _cart_coords[1] != 0; 
}

bool Mesh::has_right_neighbour() const {
  // We have a right neighbour if we are not on the last column
  return _cart_coords[1] != _dim_nodes[1] - 1; 
}
