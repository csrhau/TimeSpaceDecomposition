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


 private:
  const ConfigFile * const  _config;
  bool _debug;
  int _n_dimensions;
  std::vector<double> _min_coords;
  std::vector<double> _max_coords; 
  std::vector<int> _dim_cells; // Logical size of dimensions in cells
  std::vector<int> _padded_sizes; // _dim_cells + mesh padding
  std::vector<std::vector<double> > _cell_offsets;
  double *_u0;
  double *_u1;
  int _cell_count;
    


  void init();
};

#endif
