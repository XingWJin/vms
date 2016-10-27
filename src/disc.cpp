#include "disc.hpp"
#include "input.hpp"
#include <PCU.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <sstream>

namespace vms {

Disc::Disc(Input* in) {
  gmi_register_mesh();
  const char* g = in->geom_file.c_str();
  const char* m = in->mesh_file.c_str();
  mesh = apf::loadMdsMesh(g,m);
  apf::verify(mesh);
  dim = mesh->getDimension();
  apf::Numbering* n = apf::numberOwnedNodes(mesh,"n");
  num_owned_nodes = apf::countNodes(n);
  num_total_nodes = PCU_Add_Long(num_owned_nodes);
  numbering = apf::makeGlobal(n);
  out_file = in->out_file;
}

Disc::~Disc() {
  if (numbering)
    apf::destroyGlobalNumbering(numbering);
  if (mesh) {
    mesh->destroyNative();
    apf::destroyMesh(mesh);
  }
}

void Disc::write(const int i) {
  std::ostringstream oss;
  oss << out_file << "_" << i;
  apf::writeVtkFiles(oss.str().c_str(), mesh);
}

}
