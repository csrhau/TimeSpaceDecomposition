#ifndef MESH_H
#define MESH_H

class Mesh {
 public:
  virtual ~Mesh() = 0;
  virtual void advance() = 0;
  virtual void reflect_boundary(int boundary_) = 0;
  virtual double * get_u0() = 0;
  virtual double * get_u1() = 0;
  virtual double get_del_x() const = 0;
  virtual double get_del_y() const = 0;
  // These get the coordinates of a row, column in our matrix
  // NOTE: This includes the padded/ghost cells!
  virtual double get_x_coord(int j_) const = 0; 
  virtual double get_y_coord(int i_) const = 0; 
  virtual int get_node_core_row_count() const = 0;
  virtual int get_node_core_col_count() const = 0;
  virtual int get_node_augmented_row_count() const = 0;
  virtual int get_node_augmented_col_count() const = 0;
  virtual int get_node_core_cell_count() const = 0;
  virtual int get_node_augmented_cell_count() const = 0;

};

inline Mesh::~Mesh() {}
#endif
