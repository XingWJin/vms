#ifndef vms_solver_hpp
#define vms_solver_hpp

#include "input.hpp"
#include "mesh.hpp"
#include "linear_algebra.hpp"

namespace vms {

class Solver {
  public:
    Solver(Input* in);
    ~Solver();
    Mesh* get_mesh() {return mesh;}
    void solve();
  private:
    Input* in;
    Mesh* mesh;
    LinearAlgebra* la;
};

}

#endif
