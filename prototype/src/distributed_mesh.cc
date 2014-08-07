#include "distributed_mesh.h"

#include "config_file.h"

DistributedMesh::DistributedMesh(const ConfigFile& config_,
                                 const MPI_Comm cart_comm_,
                                 const std::vector<int>& dim_nodes_) : Mesh(config_),
                                                                       _cart_comm(cart_comm_),
                                                                       _dim_nodes(dim_nodes_) {
  MPI_Comm_rank(_cart_comm, &_cart_rank);  
  _cart_coords.resize(2, 0);
  MPI_Cart_coords(_cart_comm, _cart_rank, 2, &_cart_coords[0]);
  // Compute ranks of neighbours
  _neighbour_rank_or_neg.resize(4, -1);
  int coords[2];
  if (has_top_neighbour()) {
    coords[0] = get_node_row() - 1;
    coords[1] = get_node_col();
    MPI_Cart_rank(_cart_comm, coords, &_neighbour_rank_or_neg[TOP]);
  }
  if (has_bottom_neighbour()) {
    coords[0] = get_node_row() + 1;
    coords[1] = get_node_col();
    MPI_Cart_rank(_cart_comm, coords, &_neighbour_rank_or_neg[BOTTOM]);
  }
  if (has_left_neighbour()) {
    coords[0] = get_node_row();
    coords[1] = get_node_col() - 1;
    MPI_Cart_rank(_cart_comm, coords, &_neighbour_rank_or_neg[LEFT]);
  }
  if (has_right_neighbour()) {
    coords[0] = get_node_row();
    coords[1] = get_node_col() + 1;
    MPI_Cart_rank(_cart_comm, coords, &_neighbour_rank_or_neg[RIGHT]);
  }
}

int DistributedMesh::get_neighbour_rank(int neighbour_) const {
  return _neighbour_rank_or_neg.at(neighbour_);
}

int DistributedMesh::get_node_row() const {
  return _cart_coords.at(0);
}

int DistributedMesh::get_node_col() const {
  return _cart_coords.at(1);
}

int DistributedMesh::get_vertical_nodes_count() const {
  return _dim_nodes.at(0);
}

int DistributedMesh::get_horizontal_nodes_count() const {
  return _dim_nodes.at(1);
}

bool DistributedMesh::has_top_neighbour() const {
  // We have a top neighbour if we are not on the first row
  return get_node_row() != 0; 
}

bool DistributedMesh::has_bottom_neighbour() const {
  // We have a bottom neighbour if we are not on the last row
  return get_node_row() != get_vertical_nodes_count() - 1; 
}

bool DistributedMesh::has_left_neighbour() const {
  // We have a left neighbour if we are not on the first column
  return get_node_col() != 0; 
}

bool DistributedMesh::has_right_neighbour() const {
  // We have a right neighbour if we are not on the last column
  return get_node_col() != get_horizontal_nodes_count() - 1; 
}
