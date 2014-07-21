#ifndef CALCULATION_H
#define CALCULATION_H

class ConfigFile;
class Mesh;
class Calculation {
 public: 
  Calculation(const ConfigFile * const config_,
              Mesh * const mesh_);
  ~Calculation();
  void step();
  int n_dimensions() const;
 private:
  int _n_dimensions;
  const ConfigFile * const  _config;
  Mesh * const _mesh;
  void init();
};

#endif
