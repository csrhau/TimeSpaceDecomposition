#include <mpi.h>
#include <iostream>
#include <string>
#include <stdexcept>

#include "config_file.h"
#include "driver.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  if (argc != 2) {
    std::cerr << "Usage: deqn <filename>" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  try {
    ConfigFile config(argv[1]);
    Driver driver(config);
    driver.run();
  } catch (std::logic_error& ex) {
    std::cerr << "Exception thrown: " << ex.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPI_Finalize();
  return 0;
}
