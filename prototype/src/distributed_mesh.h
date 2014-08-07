#ifndef DISTRIBUTED_MESH_H
#define DISTRIBUTED_MESH_H

#include "mesh.h"

#include <mpi.h>
#include <vector>

namespace {
  enum Boundary {
    TOP = 0,
    BOTTOM = 1,
    LEFT = 2,
    RIGHT = 3,
  };
}

class ConfigFile;
class DistributedMesh : public Mesh {
 public:
  DistributedMesh(const ConfigFile& config_, 
                  const MPI_Comm cart_comm_,
                  const std::vector<int>& dim_nodes_);
  virtual ~DistributedMesh() = 0;
  // MPI specific things
  int get_neighbour_rank(int rank_) const;
  int get_node_row() const;
  int get_node_col() const;
  int get_vertical_nodes_count() const;
  int get_horizontal_nodes_count() const;
  bool has_top_neighbour() const;
  bool has_bottom_neighbour() const;
  bool has_left_neighbour() const;
  bool has_right_neighbour() const;
 protected:
  const MPI_Comm _cart_comm;
 private:
  std::vector<int> _cart_coords;
  const std::vector<int>& _dim_nodes;
  int _cart_rank;
  std::vector<int> _neighbour_rank_or_neg;
  virtual void exchange_boundaries() = 0;
};

inline DistributedMesh::~DistributedMesh() {}
#endif
