#include "config_file.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

ConfigFile::ConfigFile(std::string& filename_) : _filename(filename_) {
  std::ifstream ifs(filename_.c_str());
  if (!ifs.good()) {
    std::cerr << "Input file " << filename_ << " not found!" << std::endl;
    exit(1);
  }
  while (!ifs.eof()) {
    std::string line;
    std::getline(ifs, line);
    std::istringstream iss(line);

    // Get key out. First text on a line - ignore blank lines and comments
    std::string key;
    iss >> key;
    if (key.empty() || key[0] == '#') {
      continue;
    }
    std::string value;
    std::getline(iss, value);
    // Remove comments
    size_t pos = value.find_first_of('#');
    if (pos != std::string::npos) {
      value.erase(pos);
    }
    // RTRIM
    pos = value.find_last_not_of(' ');
    if (pos != std::string::npos) { 
      value.erase(pos + 1);
      pos = value.find_first_not_of(' ');
      if (pos != std::string::npos) {
        value.erase(0, pos);
      }
    } else {
      value.erase(value.begin(), value.end());
    }  
    if (value.empty()) {
     std::cerr << "Invalid configuration in " << _filename << ", Key "
               << key << " was specified without a value!" << std::endl;
     exit(1);
    }
    if(_config_mapping.find(key) != _config_mapping.end()) {
            std::cerr << "Duplicate key " << key << " in input file" << std::endl;
            exit(1);
    }
    _config_mapping[key] = value;
  }
}

ConfigFile::~ConfigFile() {}
