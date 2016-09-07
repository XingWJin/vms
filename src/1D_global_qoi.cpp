#include "control.hpp"
#include "input.hpp"
#include "solver.hpp"
#include "estimator.hpp"

namespace {

double k;
double a;

double f(apf::Vector3 const& x) {
  return 1.0;
}

double q(apf::Vector3 const& x) {
  return 1.0;
}

double u(apf::Vector3 const& x) {
  return (1.0/a)*(x[0]-(exp((a/k)*(x[0]-1.0))-exp(-a/k))/(1.0-exp(-a/k)));
}

vms::Input in = {
  /* spatial dimension=*/       1,
  /* num 1D grid points=*/      0,
  /* simplical elements=*/      false,
  /* diffusive coefficient=*/   0.0,
  /* advective coefficient=*/   apf::Vector3(0.0,0.0,0.0),
  /* forcing function=*/        f,
  /* qoi function=*/            q,
  /* exact solution=*/          u,
  /* exact qoi=*/               0.4990,
  /* output name=*/             ""};

static void print_usage(char const* exe) {
  vms::print("usage:");
  vms::print("%s <num elems> <k> <a> <output name>");
}

static void check_args(int argc, char** argv) {
  if (argc != 5) {
    print_usage(argv[0]);
    vms::fail("incorrect number of arguments");
  }
}

static void setup_input(vms::Input* in, char** argv) {
  in->num_elems = atoi(argv[1]);
  in->k = atof(argv[2]);
  in->a[0] = atof(argv[3]);
  in->output_name = std::string(argv[4]);
  k = in->k;
  a = in->a[0];
  vms::print("running with the inputs:");
  vms::print(" num 1D grid elems: %d", in->num_elems);
  vms::print(" k:                 %e", in->k);
  vms::print(" a:                 %e", in->a[0]);
  vms::print(" output name:       %s", in->output_name.c_str());
}

static void run_example(vms::Input* in) {
  vms::Solver solver(in);
  solver.solve();
  vms::Estimator estimator(in, solver.get_mesh());
  estimator.estimate();
  solver.get_mesh()->write(in->output_name.c_str());
  estimator.summarize();
}

}

int main(int argc, char** argv) {
  vms::initialize();
  check_args(argc, argv);
  setup_input(&in, argv);
  run_example(&in);
  vms::finalize();
}
