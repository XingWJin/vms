#include "mesh.hpp"
#include "control.hpp"

#include <PCU.h>
#include <apfBox.h>
#include <parma.h>

namespace vms {

static apf::Mesh2* create_mesh(int d, int n, bool is) {
  int nx = n;
  int ny = 0; 
  int nz = 0;
  if (d > 1) ny = n;
  if (d > 2) nz = n;
  double wx = 1.0;
  double wy = 0.0;
  double wz = 0.0;
  if (d > 1) wy = 1.0;
  if (d > 2) wz = 1.0;
  apf::Mesh2* m = apf::makeMdsBox(nx,ny,nz,wx,wy,wz,is);
  apf::verify(m);
  return m;
}

Mesh::Mesh(int d, int n, bool is) :
  dim(d),
  num_nodes(0),
  numbering(0) {
    mesh = create_mesh(d,n,is);
    numbering = apf::makeGlobal(numberOwnedNodes(mesh,"numbering"));
    num_nodes = apf::countNodes(numbering);
}

Mesh::~Mesh() {
  if (numbering)
    apf::destroyGlobalNumbering(numbering);
  if (mesh) {
    mesh->destroyNative();
    apf::destroyMesh(mesh);
  }
}

double Mesh::get_mesh_size(apf::MeshEntity* e) {
  double h = 0.0;
  apf::Downward edges;
  int ne = mesh->getDownward(e,1,edges);
  for (unsigned i=0; i < ne; ++i)
    h += std::max(h, apf::measure(mesh,edges[i]));
  return h/ne;
}

void Mesh::write(char const* name) {
  apf::writeVtkFiles(name,mesh);
}

}
