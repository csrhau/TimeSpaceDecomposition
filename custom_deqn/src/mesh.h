#ifndef MESH_H
#define MESH_H

#include <vector>
#include <mpi.h>

class ConfigFile;
class Mesh {
 public:
  Mesh(const ConfigFile& config_,
       MPI_Comm cart_comm_, 
       const std::vector<int>& dim_nodes_);
  ~Mesh();
  void step();
  int n_dimensions() const;
  double * get_u0();
  double * get_u1();
  int get_cell_count();
  double displacement(int dimension, int cell) const;
  const std::vector<int>& get_padded_sizes() const;
  int get_from_index(int dim_) const;
  int get_to_index(int dim_) const;
  double get_dimension_delta(int dim_) const;
  int get_internal_cells(int dim_) const;
  int get_dimension_span(int dim_) const;
  void update_boundaries();
  bool has_prev_node(int dim_) const;
  bool has_next_node(int dim_) const;

 private:
  const ConfigFile& _config;
  bool _debug;
  int _n_dimensions;
  std::vector<double> _min_coords;
  std::vector<double> _max_coords;
  std::vector<double> _dimension_deltas;
  std::vector<int> _from_indexes;
  std::vector<int> _to_indexes;
  std::vector<int> _dim_cells; // Logical size of dimensions in cells
  std::vector<int> _padded_sizes; // _dim_cells + mesh padding
  std::vector<std::vector<double> > _cell_offsets;
  double *_u0;
  double *_u1;
  int _cell_count;
  void parse_config();
  void init();
  void reflect_boundary_impl(int boundary_);
  void reflect_boundary(int dim_, bool lower_boundary_);
  // MPI stuff
  MPI_Comm _cart_comm;
  const std::vector<int>& _dim_nodes;
  int _cart_rank;
  std::vector<int> _cart_coords;
};

#endif
