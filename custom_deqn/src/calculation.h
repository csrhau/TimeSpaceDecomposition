#ifndef CALCULATION_H
#define CALCULATION_H

class ConfigFile;
class Calculation {
 public: 
  Calculation(const ConfigFile& config_);
  ~Calculation();
 private:
  const ConfigFile& _config;
  void init();
  void step();
};

#endif
