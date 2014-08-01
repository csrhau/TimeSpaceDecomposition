#ifndef STATIC_BLOCKING_MESH_H
#define STATIC_BLOCKING_MESH_H
#include "distributed_mesh.h" 
#include <mpi.h>
#include <vector>
class ConfigFile;
class StaticBlockingMesh : public DistributedMesh {
 public:
  StaticBlockingMesh(const ConfigFile& config_,
       MPI_Comm cart_comm_,
       const std::vector<int>& dim_nodes_);
  virtual ~StaticBlockingMesh();

  void advance();
  void reflect_boundary(int boundary_);
  double * get_u0();
  double * get_u1();
  int get_node_core_row_count() const;
  int get_node_core_col_count() const;
  int get_node_augmented_row_count() const;
  int get_node_augmented_col_count() const;
  int get_node_core_cell_count() const;
  int get_node_augmented_cell_count() const;
  double get_world_core_row_count() const;
  double get_world_core_col_count() const;
  // x, y signify physical (as opposed to logical) coordinates
  double get_del_y() const;
  double get_del_x() const;
  double get_core_row_y(int row_) const;
  double get_core_col_x(int col_) const;
  double get_y_coord(int row_) const;
  double get_x_coord(int col_) const;
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
  // core meaning not including boundaries, ghosts
  int _world_core_row_count;
  int _world_core_col_count;
  double _world_height; // Corresponds to rows
  double _world_width;  // Corresponds to cols
  int _node_core_row_count;
  int _node_core_col_count;
  double _core_origin_x; 
  double _core_origin_y;
  // augmented meaning boundaries, ghosts are included
  int _node_augmented_row_count;
  int _node_augmented_col_count;
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
