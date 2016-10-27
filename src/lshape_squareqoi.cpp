#include "control.hpp"
#include "input.hpp"
#include "disc.hpp"
#include "solver.hpp"
#include "estimator.hpp"

namespace {

static const double k = 0.001;
static const double a0 = -1.0;
static const double a1 = 1.0;
static const double Ju = 0.0;

double f(apf::Vector3 const&) {
  return 1.0;
}

double q(apf::Vector3 const& x) {
  if ( (x[0] <= -0.5) && (x[1] >= 0.5))
    return 1.0;
  else
    return 0.0;
}

vms::Input in = {
  /* diffusive coefficient */   k,
  /* advective coefficient */   apf::Vector3(a0,a1,0.0),
  /* forcing function */        f,
  /* dual forcing function */   q,
  /* exact qoi value */         Ju,
  /* adapt method */            vms::SPR,
  /* geom file */               "",
  /* mesh file */               "",
  /* out file */                ""
};

void solve_primal(vms::Disc* disc) {
  vms::Solver primal(&in, disc, false);
  primal.solve();
}

void solve_dual(vms::Disc* disc) {
  vms::Solver dual(&in, disc, true);
  dual.solve();
}

void estimate_error(vms::Disc* disc) {
  vms::Estimator error(&in, disc);
  error.estimate();
  error.summarize();
}

void run() {
  vms::Disc disc(&in);
  solve_primal(&disc);
  solve_dual(&disc);
  estimate_error(&disc);
  disc.write(0);
}

}

int main(int argc, char** argv) {
  vms::initialize();
  in.geom_file = argv[1];
  in.mesh_file = argv[2];
  in.out_file = argv[3];
  run();
  vms::finalize();
}
