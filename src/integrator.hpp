#ifndef vms_integrator_hpp
#define vms_integrator_hpp

#include "input.hpp"
#include "disc.hpp"
#include <apf.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>

namespace vms {

class Integrator : public apf::Integrator {
  public:
    Integrator(Input* in, Disc* d, bool is_dual);
    void inElement(apf::MeshElement* me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dv);
    apf::DynamicVector Fe;
    apf::DynamicMatrix Ke;
  private:
    Disc* disc;
    function f;
    double k;
    apf::Vector3 a;
    int num_dims;
    int num_dofs;
    apf::Element* e;
};

}

#endif
