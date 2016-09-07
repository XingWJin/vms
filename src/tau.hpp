#ifndef vms_tau_hpp
#define vms_tau_hpp

#include <apf.h>

namespace vms {

class Mesh;

double get_tau(Mesh* m, apf::MeshElement* e, double k, apf::Vector3 a);

}

#endif
