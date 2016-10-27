#include "tau.hpp"
#include "disc.hpp"

namespace vms {

static double get_mesh_size(Disc* d, apf::MeshEntity* e) {
  double h = 0.0;
  apf::Downward edges;
  apf::Mesh* m = d->get_apf_mesh();
  int ne = m->getDownward(e, 1, edges);
  for (int i=0; i < ne; ++i) {
    double length = apf::measure(m, edges[i]);
    h += length*length;
  }
  return sqrt(h/ne);
}

double get_tau(
    Disc* d,
    apf::MeshEntity* e,
    double k,
    apf::Vector3 const& a) {
  double h = get_mesh_size(d, e);
  double am = a.getLength();
  double alpha = h*am/(2.0*k);
  return h/(2.0*am)*(1.0/std::tanh(alpha)-1.0/alpha);
}

}
