#ifndef vms_mesh_hpp
#define vms_mesh_hpp

#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

namespace vms {

class Mesh {
  public:
    Mesh(int d, int n, bool is);
    ~Mesh();
    int get_dim() {return dim;}
    long get_num_nodes() {return num_nodes;}
    apf::Mesh2* get_apf_mesh() {return mesh;}
    apf::GlobalNumbering* get_numbering() {return numbering;}
    double get_mesh_size(apf::MeshEntity* e);
    void write(char const* name);
  private:
    int dim;
    long num_nodes;
    apf::Mesh2* mesh;
    apf::GlobalNumbering* numbering;
};

};

#endif
