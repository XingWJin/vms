#ifndef vms_input_hpp
#define vms_input_hpp

#include <apf.h>

namespace vms {

typedef double(*function)(apf::Vector3 const& x);

struct Input {
  unsigned spatial_dim;
  unsigned num_elems;
  bool simplical_elems;
  double k;
  apf::Vector3 a;
  function forcing_function;
  function qoi_function;
  function exact_solution;
  std::string output_name;
};

}

#endif
