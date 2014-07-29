#ifndef CALCULATION_H
#define CALCULATION_H
class ConfigFile;
class Mesh;
class Calculation {
 public:
  Calculation(const ConfigFile& config_, Mesh *mesh_);
  ~Calculation();
  void step(double dt_);
 private:
  void diffuse(double dt_);
  const ConfigFile& _config;
  Mesh * const _mesh;
};
#endif
