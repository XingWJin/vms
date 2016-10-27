#ifndef vms_input_hpp
#define vms_input_hpp

#include <apf.h>

namespace vms {

enum Method {SPR, VMS1, VMS2};

typedef double(*function)(apf::Vector3 const& x);

struct Input {
  double k;
  apf::Vector3 a;
  function forcing_function;
  function dual_forcing_function;
  double exact_qoi;
  Method adapt_method;
  std::string geom_file;
  std::string mesh_file;
  std::string out_file;
};

}

#endif
