#ifndef MESH_H
#define MESH_H
class ConfigFile;
class Mesh {
 public:
  Mesh(const ConfigFile& config_);
  virtual ~Mesh() = 0;
  virtual void advance() = 0;
  virtual void reflect_boundary(int boundary_) = 0;
  virtual double * get_u0() = 0;
  virtual double * get_u1() = 0;
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
  virtual int get_current_row_offset() const = 0;
  virtual int get_current_col_offset() const = 0;
  virtual int get_previous_row_offset() const = 0;
  virtual int get_previous_col_offset() const = 0;
  double get_world_core_row_count() const;
  double get_world_core_col_count() const;
  double get_world_height() const;
  double get_world_width() const;
  double get_del_x() const;
  double get_del_y() const;
 protected:
  const ConfigFile& _config;
 private:
  // TODO Probably actually better stored in vectors.
  int _world_core_row_count;
  int _world_core_col_count;
  double _world_height; // Corresponds to rows                                     
  double _world_width;  // Corresponds to cols  
};

inline Mesh::~Mesh() {}
#endif
