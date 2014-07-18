#ifndef EXPERIMENT_H
#define EXPERIMENT_H

class ConfigFile;

class Experiment {
 public:
  Experiment(const ConfigFile confg_);
  ~Experiment();
  void run();

 private:
  const ConfigFile& _config;

}; 

#endif
