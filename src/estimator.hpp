#ifndef vms_estimator_hpp
#define vms_estimator_hpp

#include "input.hpp"
#include "mesh.hpp"

#include <apf.h>

namespace vms {

class Estimator : public apf::Integrator {
  public:
    Estimator(Input* in, Mesh* m);
    void inElement(apf::MeshElement* me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dv);
    void estimate();
    void summarize();
  private:
    Input* in;
    Mesh* mesh;
    apf::Field* fem_solution;
    apf::Field* exact_error;
    apf::Field* approx_error;
    apf::Field* effectivities;
    apf::Field* fine_scale_solution;
    apf::Element* uh_elem;
    double Juh;
    double Jeh;
    double Jeh_bound;
    double Je_elem;
    double Jeh_elem;
};

}

#endif
