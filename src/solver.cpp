#include "solver.hpp"
#include "integrator.hpp"
#include "control.hpp"
#include <gmi.h>

namespace vms {

Solver::Solver(Input* input, Disc* d, bool is) {
  in = input;
  disc = d;
  is_dual = is;
  la = new LinAlg(d->get_num_owned_nodes(), d->get_num_total_nodes());
}

Solver::~Solver() {
  if (la) delete la;
}

static void add_to_system(
    apf::DynamicVector& Fe,
    apf::DynamicMatrix& Ke,
    apf::MeshEntity* e,
    apf::GlobalNumbering* n,
    LinAlg* la) {
  apf::NewArray<long> numbers;
  int sz = int(apf::getElementNumbers(n,e,numbers));
  la->add_to_vector(sz,&numbers[0],&Fe[0]);
  la->add_to_matrix(sz,&numbers[0],&(Ke(0,0)));
}

static void assemble_volumetric(
    Input* in, Disc* disc, LinAlg* la, bool is_dual) {
  apf::Mesh* m = disc->get_apf_mesh();
  apf::GlobalNumbering* n = disc->get_numbering();
  Integrator adv_diff(in, disc, is_dual);
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = m->begin(m->getDimension());
  while ((elem = m->iterate(elems))) {
    apf::MeshElement* me = apf::createMeshElement(m,elem);
    adv_diff.process(me);
    add_to_system(adv_diff.Fe,adv_diff.Ke,elem,n,la);
    apf::destroyMeshElement(me);
  }
  m->end(elems);
  la->synchronize();
}

static void assemble_dirichlet(Disc* disc, LinAlg* la) {
  apf::Mesh* m = disc->get_apf_mesh();
  apf::GlobalNumbering* gn = disc->get_numbering();
  gmi_model* model = m->getModel();
  gmi_ent* boundary;
  gmi_iter* boundaries = gmi_begin(model,m->getDimension()-1);
  int num_boundaries = 0;
  while ((boundary = gmi_next(model,boundaries))) {
    apf::DynamicArray<apf::Node> nodes;
    apf::ModelEntity* b = reinterpret_cast<apf::ModelEntity*>(boundary);
    int tag = m->getModelTag(b);
    if ( (tag == 16) || (tag == 18) || (tag == 20) || (tag == 22) )
      continue;
    apf::getNodesOnClosure(m,b,nodes,m->getShape());
    int n_nodes = nodes.getSize();
    long rows[n_nodes];
    for (int n=0; n < n_nodes; ++n)
      rows[n] = apf::getNumber(gn, nodes[n]);
    la->diag_matrix_rows(n_nodes, rows);
    la->zero_vector_rows(n_nodes, rows);
    num_boundaries++;
  }
  gmi_end(model, boundaries);
  la->synchronize();
  print("homogeneous dbcs applied to %d boundaries", num_boundaries);
}

void Solver::solve() {
  assemble_volumetric(in,disc,la,is_dual);
  assemble_dirichlet(disc,la);
  la->solve();
  if (!is_dual)
    la->attach(disc, "uh");
  else
    la->attach(disc, "zh");
}

}
