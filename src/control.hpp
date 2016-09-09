#ifndef vms_control_hpp
#define vms_control_hpp

namespace vms {

void initialize();
void finalize();
double time();
void print(char const* format, ...);
void fail(char const* format, ...) __attribute__((noreturn));

}

#define ASSERT(c) \
    ((c) ? ((void)0) : \
     vms::fail("assertion (%s) failed at %s:%d\n",#c,__FILE__,__LINE__))

#endif
