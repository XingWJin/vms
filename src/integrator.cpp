#include "integrator.hpp"
#include "tau.hpp"

namespace vms {

Integrator::Integrator(Input* in, Disc* d, bool is_dual) :
  apf::Integrator(1) {
  disc = d;
  f = in->forcing_function;
  k = in->k;
  a = in->a;
  num_dims = disc->get_dim();
  num_dofs = 0;
  e = 0;
  if (is_dual) {
    f = in->dual_forcing_function;
    for (int i=0; i < num_dims; ++i)
      a[i] = -a[i];
  }
}

void Integrator::inElement(apf::MeshElement* me) {
  apf::Field* c = disc->get_apf_mesh()->getCoordinateField();
  e = apf::createElement(c,me);
  num_dofs = apf::countNodes(e);
  Fe.setSize(num_dofs);
  Ke.setSize(num_dofs,num_dofs);
  for (int i=0; i < num_dofs; ++i)
    Fe(i) = 0.0;
  for (int i=0; i < num_dofs; ++i)
    for (int j=0; j < num_dofs; ++j)
      Ke(i,j) = 0.0;
}

void Integrator::outElement() {
  apf::destroyElement(e);
}

void Integrator::atPoint(apf::Vector3 const& p, double w, double dv) {

  apf::Mesh* m = disc->get_apf_mesh();
  apf::MeshElement* me = apf::getMeshElement(e);
  apf::MeshEntity* ent = apf::getMeshEntity(me);

  apf::FieldShape* s = apf::getShape(m->getCoordinateField());

  apf::NewArray<double> BF;
  apf::getBF(s,me,p,BF);

  apf::NewArray<apf::Vector3> gBF;
  apf::getGradBF(s,me,p,gBF);

  apf::Vector3 x;
  apf::mapLocalToGlobal(me,p,x);

  for (int i=0; i < num_dofs; ++i) {
    Fe(i) += f(x)*BF[i];
    Fe(i) *= w*dv;
  }

  for (int i=0; i < num_dofs; ++i) {
    for (int j=0; j < num_dofs; ++j) {
      for (int d=0; d < num_dims; ++d) {
        Ke(i,j) += k*gBF[j][d]*gBF[i][d];
        Ke(i,j) += a[d]*gBF[j][d]*BF[i];
      }
      Ke(i,j) *= w*dv;
    }
  }

  double tau = get_tau(disc, ent, k, a);

  for (int i=0; i < num_dofs; ++i) {
    double adv = 0.0;
    for (int d=0; d < num_dims; ++d)
      adv += a[d]*gBF[i][d];
    Fe(i) += tau*f(x)*adv*w*dv;
  }

  for (int i=0; i < num_dofs; ++i) {
    for (int j=0; j < num_dofs; ++j) {
      double adv1 = 0.0;
      double adv2 = 0.0;
      for (int d=0; d < num_dims; ++d) {
        adv1 += a[d]*gBF[j][d];
        adv2 += a[d]*gBF[i][d];
      }
      Ke(i,j) += tau*adv1*adv2*w*dv;
    }
  }

}

}
