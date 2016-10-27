#ifndef vms_solver_hpp
#define vms_solver_hpp

#include "input.hpp"
#include "disc.hpp"
#include "linalg.hpp"

namespace vms {

class Solver {
  public:
    Solver(Input* in, Disc* d, bool is);
    ~Solver();
    void solve();
  private:
    Input* in;
    Disc* disc;
    LinAlg* la;
    bool is_dual;
};

}

#endif
