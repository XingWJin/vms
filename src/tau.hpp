#ifndef vms_tau_hpp
#define vms_tau_hpp

#include <apf.h>

namespace vms {

class Disc;

double get_tau(
    Disc* d,
    apf::MeshEntity* e,
    double k,
    apf::Vector3 const& a);

}

#endif
