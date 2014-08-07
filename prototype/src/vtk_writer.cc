#include "vtk_writer.h"

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
      ofs.close();
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
    const int core_rows = _mesh->get_node_core_row_count();
    const int core_cols = _mesh->get_node_core_col_count();
    const int core_cells = _mesh->get_node_core_cell_count();  
    const int row_offset = _mesh->get_current_row_offset();
    const int col_offset = _mesh->get_current_col_offset();

    const int horizontal_points = core_rows + 1;
    const int vertical_points = core_cols + 1;
    file << "DIMENSIONS " << horizontal_points << " " << vertical_points 
         << " 1" << std::endl;

    file << "X_COORDINATES " << horizontal_points << " float" << std::endl;
    for(int j = 0; j < horizontal_points; ++j) {
        file << _mesh->get_x_coord(j) << " ";
    }
    file << std::endl;
    file << "Y_COORDINATES " << vertical_points << " float" << std::endl;
    for(int i = 0; i < vertical_points; ++i) {
        file << _mesh->get_y_coord(i) << " ";
    }
    file << std::endl;

    file << "Z_COORDINATES 1 float" << std::endl;
    file << "0.0000" << std::endl;

    file << "CELL_DATA " << core_cells << std::endl;

    file << "FIELD FieldData 1" << std::endl;

    file << "u 1 " << core_cells << " double" <<  std::endl;

    int x_span = _mesh->get_node_augmented_col_count();
    double *u0 = _mesh->get_u0();
    // N.B. Deal with padding
    for (int i = row_offset; i < core_rows + row_offset; ++i) {
      for (int j = col_offset; j < core_cols + col_offset; ++j) {
        file << u0[i * x_span + j] << " ";
      }
      file << std::endl;
    }
    file.close();
}
