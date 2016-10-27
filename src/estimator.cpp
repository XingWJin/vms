#include "estimator.hpp"
#include "tau.hpp"
#include "control.hpp"
#include <apfShape.h>

namespace vms {

Estimator::Estimator(Input* in, Disc* d) :
  apf::Integrator(3) {
  disc = d;
  k = in->k;
  a = in->a;
  f = in->forcing_function;
  q = in->dual_forcing_function;
  apf::Mesh* m = disc->get_apf_mesh();
  uh = m->findField("uh");
  zh = m->findField("zh");
  apf::FieldShape* s = apf::getConstant(disc->get_dim());
  error1 = apf::createField(m, "Jeh1", apf::SCALAR, s);
  error2 = apf::createField(m, "Jeh2", apf::SCALAR, s);
  error1_bound = apf::createField(m, "Jeh1_bound", apf::SCALAR, s);
  error2_bound = apf::createField(m, "Jeh2_bound", apf::SCALAR, s);
  Ju = in->exact_qoi;
  Juh = 0.0;
  Jeh1 = 0.0;
  Jeh2 = 0.0;
  Jeh1_bound = 0.0;
  Jeh2_bound = 0.0;
  uh_elem = 0;
  zh_elem = 0;
}

void Estimator::inElement(apf::MeshElement* me) {
  Jeh1_elem = 0.0;
  Jeh2_elem_term1 = 0.0;
  Jeh2_elem_term2 = 0.0;
  uh_elem = apf::createElement(uh, me);
  zh_elem = apf::createElement(zh, me);
}

void Estimator::outElement() {
  double a = Jeh1_elem;
  double b = Jeh2_elem_term1 + Jeh2_elem_term2;
  double c = std::abs(Jeh1_elem);
  double d = std::abs(Jeh2_elem_term1) + std::abs(Jeh2_elem_term2);
  apf::MeshEntity* e = apf::getMeshEntity(uh_elem);
  apf::setScalar(error1, e, 0, a);
  apf::setScalar(error2, e, 0, b);
  apf::setScalar(error1_bound, e, 0, c);
  apf::setScalar(error2_bound, e, 0, d);
  Jeh1 += a;
  Jeh2 += b;
  Jeh1_bound += c;
  Jeh2_bound += d;
  apf::destroyElement(uh_elem);
  apf::destroyElement(zh_elem);
}

void Estimator::atPoint(apf::Vector3 const& p, double w, double dv) {
  apf::MeshEntity* mesh_ent = apf::getMeshEntity(uh_elem);
  apf::MeshElement* mesh_elem = apf::getMeshElement(uh_elem);
  apf::Vector3 x;
  apf::mapLocalToGlobal(mesh_elem, p, x);
  double u = apf::getScalar(uh_elem, p);
  double z = apf::getScalar(zh_elem, p);
  apf::Vector3 gu;
  apf::Vector3 gz;
  apf::getGrad(uh_elem, p, gu);
  apf::getGrad(zh_elem, p, gz);
  double tau = get_tau(disc, mesh_ent, k, a);
  double Ru = f(x) - (a*gu);
  double Rz = q(x) + (a*gz);
  double uprime = tau*Ru;
  double zprime = tau*Rz;
  double Lstarz = -(a*gz);
  Juh += q(x)*u*w*dv;
  Jeh1_elem += q(x)*uprime*w*dv;
  Jeh2_elem_term1 += zprime*Ru*w*dv;
  Jeh2_elem_term2 += Lstarz*uprime*w*dv;
}

void Estimator::estimate() {
  double t0 = time();
  this->process(disc->get_apf_mesh());
  double t1 = time();
  print("error estimated in %f seconds", t1-t0);
}

void Estimator::summarize() {
}

}
