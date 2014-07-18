#ifndef INPUT_FILE_H
#define INPUT_FILE_H

#include <string>
#include <map>

class ConfigFile {
 public:
   ConfigFile(std::string& filename_);
   ~ConfigFile();

 private:
   std::string _filename;
   std::map<std::string, std::string> _config_mapping;

};
#endif
