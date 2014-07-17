#ifndef CFD_ARGUMENTS_H
#define CFD_ARGUMENTS_H
#include <boost/program_options.hpp>
#include <string>

namespace argparse {
namespace po = boost::program_options;

struct Parameters {
  std::string config_file;
  int ts; // Number of timesteps to simulate
  bool v; // verbose
};

// Function used to ensure mandatory options are specified
void required_option(const po::variables_map& vm_,
                     const po::options_description& desc_,
                     const char* opt1_);

// Function used to prevent conflicting options being specified
void conflicting_options(const po::variables_map& vm_,
                         const po::options_description& desc_,
                         const char* opt1_, const char* opt2_);

// Function defines and parses command line arguments
void parse_args(int argc_, char *argv_[], Parameters& params_);

} // namespace argparse
#endif
