#ifndef vms_size_hpp
#define vms_size_hpp

#include "apf.h"

namespace vms {

apf::Field* get_iso_target_size(
    apf::Field* e,
    size_t t,
    std::string name = "size");

apf::Field* get_min_size(
    apf::Field* s1,
    apf::Field* s2,
    apf::Field* s3);

}

#endif
