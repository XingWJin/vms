#ifndef vms_disc_hpp
#define vms_disc_hpp

#include <apf.h>
#include <apfMesh2.h>
#include <apfNumbering.h>

namespace vms {

struct Input;

class Disc {
  public:
    Disc(Input* in);
    ~Disc();
    int get_dim() {return dim;}
    long get_num_owned_nodes() {return num_owned_nodes;}
    long get_num_total_nodes() {return num_total_nodes;}
    apf::Mesh2* get_apf_mesh() {return mesh;}
    apf::GlobalNumbering* get_numbering() {return numbering;}
    void update();
    void write(const int i);
    void write_pvd(const int end);
  private:
    int dim;
    long num_owned_nodes;
    long num_total_nodes;
    apf::Mesh2* mesh;
    apf::GlobalNumbering* numbering;
    std::string out_file;
};

}

#endif
