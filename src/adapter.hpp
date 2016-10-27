#ifndef vms_adapter_hpp
#define vms_adapter_hpp

#include <cstdlib>

namespace vms {

struct Input;
class Disc;

class Adapter {
  public:
    Adapter(Input* in, Disc* d);
    void adapt(size_t t, int i);
  private:
    Disc* disc;
};

}

#endif
