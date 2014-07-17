#include "arguments.h"
namespace argparse {

void required_option(const po::variables_map& vm_,
                     const po::options_description& desc_,
                     const char* opt1_)
{
    if (!vm_.count(opt1_)) {
      std::stringstream mesg;
      mesg << "Error parsing arguments. Mandatory argument '" << opt1_
           << "' is missing.\n\n" << desc_ << std::endl;
        throw std::logic_error(mesg.str());
    }
}

// Function used to prevent conflicting options being specified
void conflicting_options(const po::variables_map& vm_,
                         const po::options_description& desc_,
                         const char* opt1_, const char* opt2_)
{
    if (vm_.count(opt1_) && !vm_[opt1_].defaulted()
     && vm_.count(opt2_) && !vm_[opt2_].defaulted()) {
      std::stringstream mesg;
      mesg << "Error parsing arguments. Conflicting options '" << opt1_
           << "' and '" << opt2_ << "' present.\n\n" << desc_ << std::endl;
      throw std::logic_error(mesg.str());
    }
}

void parse_args(int argc_, char *argv_[], Parameters& params_) {
  po::options_description desc("Allowed Options");
  desc.add_options()("timesteps",
                      po::value<int>(&params_.ts)->default_value(1000),
                     "Simulation timesteps")
                    ("verbose",
                     po::value<bool>(&params_.v)->default_value(true),
                     "Verbose")
                    ("config",
                     po::value<std::string>(&params_.config_file),
                     "Path to configuration file");

  po::variables_map vm;
  po::store(po::parse_command_line(argc_, argv_, desc), vm);
  po::notify(vm);
  required_option(vm , desc, "config");
}

} // namespace argparse
