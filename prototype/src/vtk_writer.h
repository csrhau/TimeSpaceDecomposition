#ifndef VTKWRITER_H
#define VTKWRITER_H

#include <string>

class Mesh;
class VtkWriter {
    public:
        VtkWriter(std::string basename, Mesh* mesh_, int world_rank_, int world_size_);
        void write(int step, double time);

    private:
        std::string dump_basename;
        std::string vtk_header;
        Mesh* _mesh;
        int _world_rank;
        int _world_size;
        void writeVtk(int step, double time);
};
#endif
