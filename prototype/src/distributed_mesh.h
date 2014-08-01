#ifndef DISTRIBUTED_MESH_H
#define DISTRIBUTED_MESH_H

#include "mesh.h"

class DistributedMesh : public Mesh {
 public:
  virtual ~DistributedMesh() = 0;
  // MPI specific things
  bool has_top_neighbour() const;
  bool has_bottom_neighbour() const;
  bool has_left_neighbour() const;
  bool has_right_neighbour() const;
};

inline DistributedMesh::~DistributedMesh() {}
#endif
