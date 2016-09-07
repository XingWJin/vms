#include "estimator.hpp"
#include "tau.hpp"
#include "control.hpp"

#include <apfShape.h>

namespace vms {

Estimator::Estimator(Input* input, Mesh* m) :
  apf::Integrator(1),
  in(input),
  mesh(m),
  fem_solution(0),
  exact_error(0),
  approx_error(0),
  effectivities(0),
  fine_scale_solution(0),
  uh_elem(0),
  Juh(0.0),
  Jeh(0.0),
  Je_elem(0.0),
  Jeh_elem(0.0) {
    apf::Mesh2* msh = mesh->get_apf_mesh();
    fem_solution = msh->findField("uh");
    ASSERT(fem_solution);
    apf::FieldShape* s = apf::getConstant(mesh->get_dim());
    exact_error = apf::createField(msh, "Je", apf::SCALAR, s);
    approx_error = apf::createField(msh, "Jeh", apf::SCALAR, s);
    effectivities = apf::createField(msh, "I", apf::SCALAR, s);
    fine_scale_solution = apf::createField(msh, "uprime", apf::SCALAR, s);
}

void Estimator::inElement(apf::MeshElement* me) {
  Je_elem = 0.0;
  Jeh_elem = 0.0;
  uh_elem = apf::createElement(fem_solution,me);
}

static double compute_effectivity(double Je_elem, double Jeh_elem) {
  double I = Jeh_elem;
  if (Je_elem < 1.0e-13)
    I = 0.0;
  else
    I /= Je_elem;
  return I;
}

void Estimator::outElement() {
  apf::MeshEntity* mesh_ent = apf::getMeshEntity(uh_elem);
  apf::setScalar(exact_error, mesh_ent, 0, Je_elem);
  apf::setScalar(approx_error, mesh_ent, 0, Jeh_elem);
  double I = compute_effectivity(Je_elem, Jeh_elem);
  apf::setScalar(effectivities, mesh_ent, 0, I);
  apf::destroyElement(uh_elem);
}

void Estimator::atPoint(apf::Vector3 const& p, double w, double dv) {

  apf::MeshEntity* mesh_ent = apf::getMeshEntity(uh_elem);
  apf::MeshElement* mesh_elem = apf::getMeshElement(uh_elem);

  apf::Vector3 x;
  apf::mapLocalToGlobal(mesh_elem,p,x);

  double k = in->k;
  double f = in->forcing_function(x);
  double q = in->qoi_function(x);
  apf::Vector3 a = in->a;

  double u = in->exact_solution(x);
  double uh = apf::getScalar(uh_elem, p);
  apf::Vector3 uh_grad;
  apf::getGrad(uh_elem, p, uh_grad);
  double tau = get_tau(mesh, mesh_ent, k, a);
  double u_prime = tau*(f-a*uh_grad);
  apf::setScalar(fine_scale_solution, mesh_ent, 0, u_prime);

  Juh += q*uh*w*dv;
  Jeh += q*u_prime*w*dv;
  Jeh_bound += std::abs(q*u_prime*w*dv);
  Je_elem += q*(u-uh)*w*dv;
  Jeh_elem += q*u_prime*w*dv;
}

void Estimator::estimate() {
  double t0 = time();
  this->process(mesh->get_apf_mesh());
  double t1 = time();
  print("error estimated in %f seconds", t1-t0);
}

void Estimator::summarize() {
  double n = in->num_elems;
  double h = 1.0/n;
  double am = in->a.getLength();
  double k = in->k;
  double Pe = am/k;
  double Peh = h*am/(2.0*k);
  double Ju = in->exact_qoi;
  print("Summary");
  print("-------");
  print(" dofs =    %ld", mesh->get_num_nodes());
  print(" h =       %.5e", h);
  print(" Pe =      %.5e", Pe);
  print(" Peh =     %.5e", Peh);
  print(" J(u) =    %.5e", Ju);
  print(" J(uh) =   %.5e", Juh);
  print(" J(e) =    %.5e", Ju-Juh);
  print(" J(eh) ~   %.5e", Jeh);
  print(" |J(eh)| < %.5e", Jeh_bound);
  print(" I =       %.5e", Jeh/(Ju-Juh));
}

}
