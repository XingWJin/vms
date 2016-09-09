#include "control.hpp"
#include "input.hpp"
#include "solver.hpp"
#include "estimator.hpp"

namespace {

double k;
double a0;
double a1;
const double pi = 3.141592653589793;

double f(apf::Vector3 const& p) {
  double x = p[0];
  double y = p[1];
  double b0 = sin(pi*x);
  double b1 = cos(pi*x);
  double c0 = sin(pi*y);
  double c1 = cos(pi*y);
  double ux = pi*b1*c0;
  double uxx = -pi*pi*b0*c0;
  double uy = pi*b0*c1;
  double uyy = -pi*pi*b0*c0;
  return -k*(uxx+uyy) + a0*ux + a1*uy;
}

double q(apf::Vector3 const& x) {
  return 1.0;
}

double u(apf::Vector3 const& p) {
  double x = p[0];
  double y = p[1];
  return sin(pi*x)*sin(pi*y);
}

vms::Input in = {
  /* spatial dimension=*/       2,
  /* num 1D grid points=*/      0,
  /* simplical elements=*/      false,
  /* diffusive coefficient=*/   0.0,
  /* advective coefficient=*/   apf::Vector3(0.0,0.0,0.0),
  /* forcing function=*/        f,
  /* qoi function=*/            q,
  /* exact solution=*/          u,
  /* exact qoi=*/               4.0/(pi*pi),
  /* output name=*/             ""};

static void print_usage(char const* exe) {
  vms::print("usage:");
  vms::print("%s <num elems> <k> <a0> <a1> <output name>", exe);
}

static void check_args(int argc, char** argv) {
  if (argc != 6) {
    print_usage(argv[0]);
    vms::fail("incorrect number of arguments");
  }
}

static void setup_input(vms::Input* in, char** argv) {
  in->num_elems = atoi(argv[1]);
  in->k = atof(argv[2]);
  in->a[0] = atof(argv[3]);
  in->a[1] = atof(argv[4]);
  in->output_name = std::string(argv[5]);
  k = in->k;
  a0 = in->a[0];
  a1 = in->a[1];
  vms::print("running with the inputs:");
  vms::print(" num 1D grid elems: %d", in->num_elems);
  vms::print(" k:                 %e", in->k);
  vms::print(" a0:                %e", in->a[0]);
  vms::print(" a1:                %e", in->a[1]);
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
