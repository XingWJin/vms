#include "control.hpp"
#include "input.hpp"
#include "disc.hpp"
#include "solver.hpp"
#include "estimator.hpp"
#include "adapter.hpp"

namespace {

static const double k = 0.001;
static const double a0 = -1.0;
static const double a1 = 1.0;
static const double Ju = 1.65361;

double f(apf::Vector3 const&) {
  return 1.0;
}

double q(apf::Vector3 const&) {
  return 1.0;
}

vms::Input in = {
  /* diffusive coefficient */   k,
  /* advective coefficient */   apf::Vector3(a0,a1,0.0),
  /* forcing function */        f,
  /* dual forcing function */   q,
  /* exact qoi value */         Ju,
  /* adapt method */            vms::VMS2,
  /* geom file */               "",
  /* mesh file */               "",
  /* out file */                ""
};

static void solve_primal(vms::Disc* disc) {
  vms::Solver primal(&in, disc, false);
  primal.solve();
}

static void solve_dual(vms::Disc* disc) {
  vms::Solver dual(&in, disc, true);
  dual.solve();
}

static void estimate_error(vms::Disc* disc) {
  vms::Estimator error(&in, disc);
  error.estimate();
  error.summarize();
}

static void adapt_mesh(vms::Disc* disc, int i) {
  static int j=2;
  vms::Adapter adapter(&in, disc);
  adapter.adapt(j*100, i);
  j*=2;
}

void run() {
  vms::Disc disc(&in);
  for (int i=0; i < 15; ++i) {
    solve_primal(&disc);
    solve_dual(&disc);
    estimate_error(&disc);
    adapt_mesh(&disc, i);
  }
  disc.write_pvd(15);
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
