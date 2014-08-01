#include "dynamic_mesh.h"

#include <mpi.h>

namespace {
  // TODO, rename distributed_mesh to grid_mesh, then move this there
  enum Boundary {
    TOP = 1,
    BOTTOM = 2,
    LEFT = 3,
    RIGHT = 4,
  };
}

DynamicMesh::DynamicMesh(const ConfigFile& config_,
                         MPI_Comm cart_comm_,
                         int row_nodes_,
                         int col_nodes_) : _config(config_),
                                           _cart_comm(cart_comm_),
                                           _row_nodes(row_nodes_),
                                           _col_nodes(col_nodes_),
                                           _prograde(true) {
  MPI_Comm_rank(_cart_comm, &_cart_rank);
  // Get world inner simulation domains (not including ghost/boundary cells)
  std::vector<int> logical_dimensions = _config.get_or_default("logical_dimensions",
                                                               std::vector<int>());
  _world_inner_rows = logical_dimensions.at(0);
  _world_inner_cols = logical_dimensions.at(1);
  std::vector<double> physical_dimensions = _config.get_or_default("physical_dimensions",
                                                                   std::vector<double>());
  _world_height = physical_dimensions.at(0);
  _world_width = physical_dimensions.at(1);
  _cart_coords.resize(2, 0);
  MPI_Cart_coords(_cart_comm, _cart_rank, 2, &_cart_coords[0]);
  int elements = get_node_outer_cell_count();
  _u0 = new double[elements];
  _u1 = new double[elements];
}

DynamicMesh::~DynamicMesh() {
  delete[] _u0;
  delete[] _u1;
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

void reflect_boundary(int boundary_);
double * DynamicMesh::get_u0() { return _u0; }
double * DynamicMesh::get_u1() { return _u1; }
double DynamicMesh::get_del_y() const { return _world_height / _world_inner_rows; }
double DynamicMesh::get_del_x() const { return _world_width / _world_inner_cols; }
double get_inner_origin_y() const {
  return get_del_y() * calculate_local_offset(_cart_coords.at(0),
                                              _dim_nodes.at(0),
                                              logical_dimensions.at(0),
                                              _prograde);
}
double get_inner_origin_x() const {
  return get_del_x() * calculate_local_offset(_cart_coords.at(1),
                                              _dim_nodes.at(1),
                                              logical_dimensions.at(1),
                                              _prograde);
}

double get_outer_col_x(int col_) const {
  return get_inner_origin_x() + (col_ - 1) * get_del_x();
}

double get_outer_row_y(int row_) const {
  return get_inner_origin_y() + (row_ - 1) * get_del_y();
}

int get_node_inner_rows() const;
int get_node_inner_cols() const;
int get_node_outer_rows() const;
int get_node_outer_cols() const;
int get_node_inner_cell_count() const; // The size of the domain we are currently responsible for

int get_node_outer_cell_count() const {
  return get_node_augmented_row_count() * get_node_augmented_col_count();
};

bool DynamicMesh::has_top_neighbour() const {
  // We have a top neighbour if we are not on the first row
  return _cart_coords[0] != 0;
}

bool DynamicMesh::has_bottom_neighbour() const {
  // We have a bottom neighbour if we are not on the last row
  return _cart_coords[0] != _dim_nodes[0] - 1;
}

bool DynamicMesh::has_left_neighbour() const {
  // We have a left neighbour if we are not on the first column
  return _cart_coords[1] != 0;
}

bool DynamicMesh::has_right_neighbour() const {
  // We have a right neighbour if we are not on the last column
  return _cart_coords[1] != _dim_nodes[1] - 1;
}

int DynamicMesh::get_node_augmented_row_count() const {
  const int border_rows = 2;
  const int unaugmented_rows = calculate_local_span(_cart_coords.at(0),
                                                    _dim_nodes.at(0),
                                                    logical_dimensions.at(0));
  const int additional_rows = has_top_neighbour() ? 1 : 0;
}

int DynamicMesh::get_node_augmented_col_count() const {
  const int border_columns = 2;
  const int unaugmented_cols = calculate_local_span(_cart_coords.at(1),
                                                    _dim_nodes.at(1),
                                                    logical_dimensions.at(1));
  const int additional_cols = has_left_neighbour() ? 1 : 0;
}
