#ifndef vms_estimator_hpp
#define vms_estimator_hpp

#include "input.hpp"
#include "disc.hpp"
#include <apf.h>

namespace vms {

class Estimator : public apf::Integrator {
  public:
    Estimator(Input* in, Disc* d);
    void inElement(apf::MeshElement* me);
    void outElement();
    void atPoint(apf::Vector3 const& p, double w, double dv);
    void estimate();
    void summarize();
  private:
    Disc* disc;
    double k;
    apf::Vector3 a;
    function f;
    function q;
    apf::Field* uh;
    apf::Field* zh;
    apf::Field* error1;
    apf::Field* error1_bound;
    apf::Field* error2;
    apf::Field* error2_bound;
    apf::Element* uh_elem;
    apf::Element* zh_elem;
    double Ju;
    double Juh;
    double Jeh1;
    double Jeh2;
    double Jeh1_bound;
    double Jeh2_bound;
    double Jeh1_elem;
    double Jeh2_elem_term1;
    double Jeh2_elem_term2;
};

}

#endif
