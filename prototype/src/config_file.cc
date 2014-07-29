#include "config_file.h"

#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <fstream>


typedef std::map<std::string, std::string>::const_iterator config_iterator;

ConfigFile::ConfigFile(const char *filename_) : _filename(filename_) {
  std::ifstream ifs(filename_);
  if (!ifs.good()) {
    std::stringstream msg;
    msg << "Input File " << filename_ << " not found";
    throw std::logic_error(msg.str());
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
    std::size_t pos = value.find_first_of('#');
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
      std::stringstream msg;
      msg << "Invalid configuration in " << _filename << ", Key " << key 
          << " specified without value";
      throw std::logic_error(msg.str());
    }
    if(_config_mapping.find(key) != _config_mapping.end()) {
      std::stringstream msg;
      msg << "Duplicate key: " << key;
      throw std::logic_error(msg.str());
    }
    _config_mapping[key] = value;
  }
}

ConfigFile::~ConfigFile() {}

void ConfigFile::print_config() const {
  std::cout << "Run Config:";
  config_iterator it = _config_mapping.begin();
  config_iterator itEnd = _config_mapping.end();
  for (; it != itEnd; ++it) {
    std::cout << "\n\t" << it->first << " " << it->second;
  }
  std::cout << std::endl;
}

const char* ConfigFile::get_filename() const {
  return _filename;
}

// getter template function defined in header file.
// template <typename T>
// T ConfigFile::get_or_default(const std::string& name, const T& dfault) const;
