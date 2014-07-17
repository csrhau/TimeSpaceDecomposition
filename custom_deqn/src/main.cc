#include <iostream>

#include "arguments.h"
#include "InputFile.h"


int main(int argc, char *argv[]) {
  // Process user command line arguments
  argparse::Parameters params;
  try {
    parse_args(argc, argv, params);
  } catch (std::logic_error& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }

  if (params.v) {
    std::cout << "Simple 2D Heat Equation Solver using an explicit scheme\n";
    std::cout << "Timesteps: " << params.ts << "\n";
    std::cout << "Verbose: " << params.v << std::endl;
  }

  InputFile input(params.config_file.c_str());



  return 0;
}
