#ifndef DYNAMIC_MESH_H
#define DYNAMIC_MESH_H
#include "distributed_mesh.h" 
#include <mpi.h>
#include <vector>
class ConfigFile;
class DynamicMesh : public DistributedMesh {
 public:
  DynamicMesh(const ConfigFile& config_,
       MPI_Comm cart_comm_,
       int row_nodes_,
       int col_nodes_);

  virtual ~DynamicMesh();
  void advance();
  void reflect_boundary(int boundary_);
  double * get_u0();
  double * get_u1();
  // x, y signify physical (as opposed to logical) coordinates
  double get_del_x() const;
  double get_del_y() const;
  double get_inner_origin_x() const;
  double get_inner_origin_y() const;
  double get_outer_row_y(int row_) const;
  double get_outer_col_x(int col_) const;
  int get_node_inner_rows() const;
  int get_node_inner_cols() const;
  int get_node_outer_rows() const;
  int get_node_outer_cols() const;
  int get_node_inner_cell_count() const;
  int get_node_outer_cell_count() const;
  bool has_top_neighbour() const;
  bool has_bottom_neighbour() const;
  bool has_left_neighbour() const;
  bool has_right_neighbour() const;
  // MPI specific things
  bool has_top_neighbour() const;
  bool has_bottom_neighbour() const;
  bool has_left_neighbour() const;
  bool has_right_neighbour() const;

 private:
  // TODO rip out any unused members!
  const ConfigFile& _config;
  const MPI_Comm _cart_comm;
  int row_nodes_;
  int col_nodes_;
  bool _prograde;
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
  double *_u0;
  double *_u1;
  int _top_rank_or_neg;
  int _bottom_rank_or_neg;
  int _left_rank_or_neg;
  int _right_rank_or_neg;
  void exchange_boundaries();
  void get_node_augmented_row_count() const;
  void get_node_augmented_col_count() const;





};
#endif
