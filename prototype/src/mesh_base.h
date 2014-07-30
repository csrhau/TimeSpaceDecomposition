#ifndef MESH_BASE_H
#define MESH_BASE_H

class MeshBase {
 public:
  virtual ~MeshBase() = 0;
  virtual void advance() = 0;
  virtual void reflect_boundary(int boundary_) = 0;
  virtual double * get_u0() = 0;
  virtual double * get_u1() = 0;
  virtual double get_del_x() const = 0;
  virtual double get_del_y() const = 0;
  virtual double get_outer_col_x(int i_) const = 0;
  virtual double get_outer_row_y(int j_) const = 0;
  virtual int get_node_inner_rows() const = 0;
  virtual int get_node_inner_cols() const = 0;
  virtual int get_node_outer_rows() const = 0;
  virtual int get_node_outer_cols() const = 0;
  virtual int get_node_inner_cell_count() const = 0;
};

inline MeshBase::~MeshBase() {}
#endif
