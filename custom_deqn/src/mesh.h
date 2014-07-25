#ifndef MESH_H
#define MESH_H

#include <vector>

class ConfigFile;
class Mesh {
 public:
  Mesh(const ConfigFile * const config_);
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

 private:
  const ConfigFile * const  _config;
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
  void init();
  void reflect_boundary(int boundary_);
};

#endif
