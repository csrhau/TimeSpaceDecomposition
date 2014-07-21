#include "calculation.h"

#include "config_file.h"
#include "mesh.h"

Calculation::Calculation(const ConfigFile * const config_,
                         Mesh * const mesh_) : _config(config_), _mesh(mesh_) {
 _n_dimensions = _config->get_or_default("n_dimensions", 2);
}
Calculation::~Calculation() {}

void Calculation::init() {}
void Calculation::step() {} // Should be virtual
int  Calculation::n_dimensions() const { 
  return _n_dimensions;
}
