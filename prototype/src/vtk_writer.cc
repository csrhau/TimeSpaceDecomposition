#include "vtk_writer.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "mesh.h"

VtkWriter::VtkWriter(std::string basename, Mesh* mesh_, int world_rank_, int world_size_) :
    dump_basename(basename),
    vtk_header("# vtk DataFile Version 3.0\nvtk output\nASCII\n"),
    _mesh(mesh_), 
    _world_rank(world_rank_), 
    _world_size(world_size_) {

    if (_world_rank == 0) {
      std::ofstream ofs;
      std::stringstream fname;
      fname << dump_basename << ".visit";
      std::string file_name = fname.str();
      ofs.open(file_name.c_str());
      ofs << "!NBLOCKS "  << _world_size << std::endl;
      ofs.close();
    }
}

void VtkWriter::write(int step, double time)
{
    std::cout << "step " << step << " rank " << _world_rank << std::endl;
    // Master process writes out the .visit file to coordinate the .vtk files
    if (_world_rank == 0) {
      std::ofstream ofs;
      std::stringstream fname;
      fname << dump_basename << ".visit";

      std::string file_name = fname.str();

      // Open file in append mode
      ofs.open(file_name.c_str(), std::ofstream::out | std::ofstream::app);

      for (int rank = 0; rank < _world_size; ++rank) {
        ofs << dump_basename
             << "."
             << step
             << "."
             << rank
             << ".vtk" << std::endl;
      }
      ofs.flush();
      ofs.close();
      std::cout << "Written to file!" << std::endl;
    }
    writeVtk(step, time);
}

void VtkWriter::writeVtk(int step, double time)
{
    std::ofstream file;
    std::stringstream fname;
    fname << dump_basename
        << "."
        << step
        << "."
        << _world_rank
        << ".vtk";

    std::string file_name = fname.str();

    file.open(file_name.c_str());

    file.setf(std::ios::fixed, std::ios::floatfield);
    file.precision(8);

    file << vtk_header;

    file << "DATASET RECTILINEAR_GRID" << std::endl;
    file << "FIELD FieldData 2" << std::endl;
    file << "TIME 1 1 double" << std::endl;
    file << time << std::endl;
    file << "CYCLE 1 1 int" << std::endl;
    file << step << std::endl;
    // 2D - Note we need internals +1 due to fencepost effec
    int horizontal_points = _mesh->get_node_inner_cols() + 1;
    int vertical_points = _mesh->get_node_inner_rows() + 1;
    file << "DIMENSIONS " << horizontal_points << " " << vertical_points 
         << " 1" << std::endl;

    file << "X_COORDINATES " << horizontal_points << " float" << std::endl;
    for(int j = 0; j < horizontal_points; ++j) {
        file << _mesh->get_outer_col_x(j) << " ";
    }
    file << std::endl;
    file << "Y_COORDINATES " << vertical_points << " float" << std::endl;
    for(int i = 0; i < vertical_points; ++i) {
        file << _mesh->get_outer_row_y(i) << " ";
    }
    file << std::endl;

    file << "Z_COORDINATES 1 float" << std::endl;
    file << "0.0000" << std::endl;


    file << "CELL_DATA " << _mesh->get_node_inner_cell_count() << std::endl;

    file << "FIELD FieldData 1" << std::endl;

    file << "u 1 " << _mesh->get_node_inner_cell_count() << " double" <<  std::endl;

    int x_span = _mesh->get_node_outer_cols();
    double *u0 = _mesh->get_u0();
    // N.B. Deal with padding
    for (int i = 1; i < _mesh->get_node_inner_rows() + 1; ++i) {
      for (int j = 1; j < _mesh->get_node_inner_cols() + 1; ++j) {
        file << u0[i * x_span + j] << " ";
      }
      file << std::endl;
    }

    file.close();
}
