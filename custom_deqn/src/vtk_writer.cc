#include "vtk_writer.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "mesh.h"
#include "tools.h"

VtkWriter::VtkWriter(std::string basename, Mesh* mesh) :
    dump_basename(basename),
    vtk_header("# vtk DataFile Version 3.0\nvtk output\nASCII\n"),
    mesh(mesh)
{
    std::ofstream file;
    std::stringstream fname;

    fname << dump_basename << ".visit";

    std::string file_name = fname.str();

    file.open(file_name.c_str());

    file << "!NBLOCKS "
        << 1 << std::endl;
}

void VtkWriter::write(int step, double time)
{
    // Master process writes out the .visit file to coordinate the .vtk files
    std::ofstream file;
    std::stringstream fname;

    fname << dump_basename << ".visit";

    std::string file_name = fname.str();

    // Open file in append mode
    file.open(file_name.c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::app);

    file << dump_basename
        << "."
        << step
        << "."
        << 1
        << ".vtk" << std::endl;

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
        << 1
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
    if (mesh->n_dimensions() == 2) {
        file << "DIMENSIONS " << mesh->get_internal_cells(DIM_X)+1
            << " " << mesh->get_internal_cells(1)+1
            << " 1" << std::endl;
    } else if (mesh->n_dimensions() == 3) {
        file << "DIMENSIONS " << mesh->get_internal_cells(DIM_X)+1
            << " " << mesh->get_internal_cells(1)+1
            << " " << mesh->get_internal_cells(2)+1 << std::endl;
    }

    file << "X_COORDINATES " << mesh->get_internal_cells(DIM_X)+1 << " float" << std::endl;
    for(int i = 0; i < mesh->get_internal_cells(DIM_X)+1; ++i) {
        file << i * mesh->get_dimension_delta(DIM_X) << " ";
    }

    file << std::endl;

    file << "Y_COORDINATES " << mesh->get_internal_cells(1)+1 << " float" << std::endl;
    for(int j = 0; j < mesh->get_internal_cells(1)+1; ++j) {
        file << j * mesh->get_dimension_delta(1) << " ";
    }

    file << std::endl;

    file << "Z_COORDINATES 1 float" << std::endl;
    file << "0.0000" << std::endl;

    file << "CELL_DATA " << ((mesh->get_to_index(DIM_X) - mesh->get_from_index(DIM_X)) * (mesh->get_to_index(DIM_Y) - mesh->get_from_index(DIM_Y))) << std::endl;

    file << "FIELD FieldData 1" << std::endl;

    file << "u 1 " << ((mesh->get_to_index(DIM_X) - mesh->get_from_index(DIM_X)) * (mesh->get_to_index(DIM_Y) - mesh->get_from_index(DIM_Y))) << " double" <<  std::endl;

    int x_span = mesh->get_dimension_span(DIM_X);
    for(int j=mesh->get_from_index(DIM_Y); j < mesh->get_to_index(DIM_Y); ++j) {
        for (int i = mesh->get_from_index(DIM_X); i < mesh->get_to_index(DIM_X); ++i) {
            file << mesh->get_u0()[POLY2(i, j, x_span)] << " ";
        }
        file << std::endl;
    }

    file.close();
}
