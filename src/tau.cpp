#include "tau.hpp"
#include "mesh.hpp"

namespace vms {

double get_tau(
    Mesh* m, apf::MeshEntity* e, double k, apf::Vector3 const& a) {
  double h = m->get_mesh_size(e);
  double am = a.getLength();
  double alpha = h*am/(2.0*k);
  return h/(2.0*am)*(1.0/std::tanh(alpha)-1.0/alpha);
}

}
