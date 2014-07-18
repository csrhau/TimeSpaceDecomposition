#include "calculation.h"

#include "config_file.h"

Calculation::Calculation(const ConfigFile& config_) : _config(config_) {}
Calculation::~Calculation() {}

void Calculation::init() {}
void Calculation::step() {} // Should be virtual
