#include "disc.hpp"
#include "input.hpp"
#include <PCU.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <fstream>
#include <sstream>

namespace vms {

Disc::Disc(Input* in) {
  gmi_register_mesh();
  const char* g = in->geom_file.c_str();
  const char* m = in->mesh_file.c_str();
  mesh = apf::loadMdsMesh(g,m);
  apf::verify(mesh);
  dim = mesh->getDimension();
  numbering = 0;
  out_file = in->out_file;
  update();
}

Disc::~Disc() {
  if (numbering)
    apf::destroyGlobalNumbering(numbering);
  if (mesh) {
    mesh->destroyNative();
    apf::destroyMesh(mesh);
  }
}

void Disc::update() {
  if (numbering) apf::destroyGlobalNumbering(numbering);
  apf::Numbering* n = apf::numberOwnedNodes(mesh, "n");
  num_owned_nodes = apf::countNodes(n);
  num_total_nodes = PCU_Add_Long(num_owned_nodes);
  numbering = apf::makeGlobal(n);
  apf::synchronize(numbering);
}

void Disc::write(const int i) {
  std::ostringstream oss;
  oss << out_file << "_" << i;
  std::string name = oss.str();
  apf::writeVtkFiles(name.c_str(), mesh);
}

void Disc::write_pvd(const int end) {
  std::string pvd = out_file + ".pvd";
  std::fstream pvdf;
  pvdf.open(pvd.c_str(), std::ios::out);
  pvdf << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
  pvdf << "  <Collection>" << std::endl;
  for (int t=0; t < end; ++t) {
    std::ostringstream oss;
    oss << out_file << "_" << t;
    std::string vtu = oss.str();
    pvdf << "    <DataSet timestep=\"" << t << "\" group=\"\" ";
    pvdf << "part=\"0\" file=\"" << vtu << "/" << vtu;
    pvdf << ".pvtu\"/>" << std::endl;
  }
  pvdf << "  </Collection>" << std::endl;
  pvdf << "</VTKFile>" << std::endl;
}

}
