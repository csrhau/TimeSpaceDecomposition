#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <string>

class Mesh;
class VtkWriter {
    private:
        std::string dump_basename;

        std::string vtk_header;

        Mesh* mesh;

        void writeVtk(int step, double time);
    public:
        VtkWriter(std::string basename, Mesh* mesh);

        void write(int step, double time);
};
#endif
