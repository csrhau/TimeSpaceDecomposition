#ifndef STATIC_MESH_H
#define STATIC_MESH_H
#include "mesh.h" 
#include <mpi.h>
#include <vector>
class ConfigFile;
class StaticMesh : public Mesh {
 public:
  StaticMesh(const ConfigFile& config_,
       MPI_Comm cart_comm_,
       const std::vector<int>& dim_nodes_);
  virtual ~StaticMesh();

  void advance();
  void reflect_boundary(int boundary_);
  double * get_u0();
  double * get_u1();
  int get_node_inner_rows() const;
  int get_node_inner_cols() const;
  int get_node_outer_rows() const;
  int get_node_outer_cols() const;
  int get_node_inner_cell_count() const;
  int get_node_outer_cell_count() const;
  double get_world_inner_rows() const;
  double get_world_inner_cols() const;
  // x, y signify physical (as opposed to logical) coordinates
  double get_del_y() const;
  double get_del_x() const;
  double get_inner_row_y(int row_) const;
  double get_inner_col_x(int col_) const;
  double get_outer_row_y(int row_) const;
  double get_outer_col_x(int col_) const;
  // MPI specific things
  bool has_top_neighbour() const;
  bool has_bottom_neighbour() const;
  bool has_left_neighbour() const;
  bool has_right_neighbour() const;

 private:
  const ConfigFile& _config;
  const MPI_Comm _cart_comm;
  MPI_Datatype _col_type;
  std::vector<int> _cart_coords;
  // inner meaning not including boundaries, ghosts
  int _world_inner_rows;
  int _world_inner_cols;
  double _world_height; // Corresponds to rows
  double _world_width;  // Corresponds to cols
  int _node_inner_rows;
  int _node_inner_cols;
  double _inner_origin_x; 
  double _inner_origin_y;
  // Outer meaning boundaries, ghosts are included
  int _node_outer_rows;
  int _node_outer_cols;
  int _cart_rank;
  const std::vector<int>& _dim_nodes;
  double *_u0;
  double *_u1;
  int _top_rank_or_neg;
  int _bottom_rank_or_neg;
  int _left_rank_or_neg;
  int _right_rank_or_neg;
  void exchange_boundaries();
};
#endif
