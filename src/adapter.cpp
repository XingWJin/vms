#include "adapter.hpp"
#include "disc.hpp"
#include <ma.h>
#include <PCU.h>
#include <apf.h>
#include <spr.h>

namespace vms {

Adapter::Adapter(Input* in, Disc* d) {
  disc = d;
}

static void destroy_fields(Disc* d) {
  apf::Mesh* m = d->get_apf_mesh();
  apf::destroyField(m->findField("uh"));
  apf::destroyField(m->findField("zh"));
  apf::destroyField(m->findField("Jeh1"));
  apf::destroyField(m->findField("Jeh2"));
  apf::destroyField(m->findField("Jeh1_bound"));
  apf::destroyField(m->findField("Jeh2_bound"));
}

static apf::Field* get_spr_size_field(Disc* d, size_t t) {
  apf::Mesh* m = d->get_apf_mesh();
  apf::Field* uh = m->findField("uh");
  apf::Field* guh = spr::getGradIPField(uh, "guh", 1);
  apf::Field* size = spr::getTargetSPRSizeField(guh, t, 0.5, 2.0);
  apf::destroyField(guh);
  return size;
}

void Adapter::adapt(size_t t, int i) {
  apf::Mesh2* m = disc->get_apf_mesh();
  apf::Field* size = get_spr_size_field(disc, t);
  disc->write(i);
  destroy_fields(disc);
  ma::Input* in = ma::configure(m, size);
  in->maximumIterations = 3;
  ma::adapt(in);
  apf::destroyField(size);
  disc->update();
}

}
