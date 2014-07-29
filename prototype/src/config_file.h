#ifndef INPUT_FILE_H
#define INPUT_FILE_H

#include <string>
#include <map>
#include <sstream>
#include <vector>

namespace {
  typedef std::map<std::string, std::string>::const_iterator config_iterator;
}

class ConfigFile {
 public:
   ConfigFile(const char *filename_);
   ~ConfigFile();
   void print_config() const;
   const char* get_filename() const;

   // Config getters
   // General Case
   template <typename T>
   T get_or_default(const std::string& name, const T& dfault) const {
    config_iterator it = _config_mapping.find(name);
    if (it == _config_mapping.end()) {
      return dfault;
    }
    std::istringstream iss(it->second);
    T val;
    iss >> val;
    return val;
  }

  // Vectors
  template <typename T>
  const std::vector<T> get_or_default(const std::string& name,
                                    const std::vector<T>& dfault) const {
    config_iterator it = _config_mapping.find(name);
    if (it == _config_mapping.end()) {
      return dfault;
    }
    std::istringstream iss(it->second);
    std::vector<T> values;
    T val;
    while (iss >> val) {
      values.push_back(val);
    }
    return values;
  }

  // Overloads to allow 'true' to parse as true
  bool get_or_default(const std::string& name, const bool& dfault) const {
    config_iterator it = _config_mapping.find(name);
    if (it == _config_mapping.end()) {
      return dfault;
    }
    std::istringstream iss(it->second);
    bool val;
    iss >> std::boolalpha >> val;
    return val;
  }

  const std::vector<bool> get_or_default(const std::string& name,
                                    const std::vector<bool>& dfault) const {
    config_iterator it = _config_mapping.find(name);
    if (it == _config_mapping.end()) {
      return dfault;
    }
    std::istringstream iss(it->second);
    std::vector<bool> values;
    bool val;
    iss >> std::boolalpha;
    while (iss >> val) {
      values.push_back(val);
    }
    return values;
  }

 private:
   const char *_filename;
   std::map<std::string, std::string> _config_mapping;

};
#endif
