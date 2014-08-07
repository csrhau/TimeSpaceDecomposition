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
  int get_current_row_offset() const;
  int get_current_col_offset() const;
  int get_previous_row_offset() const;
  int get_previous_col_offset() const;

  // x, y signify physical (as opposed to logical) coordinates
  double get_core_row_y(int row_) const;
  double get_core_col_x(int col_) const;
  double get_y_coord(int row_) const;
  double get_x_coord(int col_) const;

 private:
  double *_u0;
  double *_u1;
  MPI_Datatype _col_type;
  // core meaning not including boundaries, ghosts
  int _world_core_col_count;
  int _node_core_row_count;
  int _node_core_col_count;
  double _core_origin_x; 
  double _core_origin_y;
  // augmented meaning boundaries, ghosts are included
  int _node_augmented_row_count;
  int _node_augmented_col_count;
  void exchange_boundaries();
};
#endif
