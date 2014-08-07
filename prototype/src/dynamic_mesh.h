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
       const std::vector<int>& dim_nodes_);

  virtual ~DynamicMesh();
  void advance();
  void reflect_boundary(int boundary_);
  double * get_u0();
  double * get_u1();
  double get_core_origin_x() const;
  double get_core_origin_y() const;
  int get_current_row_offset() const;
  int get_current_col_offset() const;
  int get_previous_row_offset() const;
  int get_previous_col_offset() const;
  double get_y_coord(int row_) const;
  double get_x_coord(int col_) const;
  int get_node_core_row_count() const;
  int get_node_core_col_count() const;
  int get_node_augmented_row_count() const;
  int get_node_augmented_col_count() const;
  int get_node_core_cell_count() const;
  int get_node_augmented_cell_count() const;

 private:
  // TODO rip out any unused members!
  bool _prograde;
  double *_u0;
  double *_u1;
  // core meaning not including boundaries, ghosts
  int _node_core_row_count;
  int _node_core_col_count;
  double _core_origin_x; // Origin at start of prograde
  double _core_origin_y; // Origin at start of prograde step
  // augmented meaning boundaries, ghosts are included
  int _node_augmented_row_count;
  int _node_augmented_col_count;
  void exchange_boundaries();
};
#endif
