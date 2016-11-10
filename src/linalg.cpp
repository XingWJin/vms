#include "linalg.hpp"
#include "disc.hpp"
#include "control.hpp"

#define CALL(function) ASSERT(0 == (function))

namespace vms {

LinAlg::LinAlg(int n, long N) {
  PetscInt NN = PetscInt(N);
  CALL(VecCreateMPI(PETSC_COMM_WORLD,n,NN,&b));
  CALL(VecSetOption(b,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE));
  CALL(MatCreateAIJ(PETSC_COMM_WORLD,n,n,NN,NN,300,PETSC_NULL,300,PETSC_NULL,&A));
  CALL(MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE));
  CALL(KSPCreate(PETSC_COMM_WORLD,&solver));
  CALL(KSPSetTolerances(solver,1.0e-10,1.0e-10,PETSC_DEFAULT,1000));
  CALL(VecDuplicate(b,&x));
}

LinAlg::~LinAlg() {
  CALL(MatDestroy(&A));
  CALL(VecDestroy(&x));
  CALL(VecDestroy(&b));
  CALL(KSPDestroy(&solver));
}

void LinAlg::add_to_vector(int sz, long* rows, double* vals) {
  PetscInt r[sz];
  PetscInt s = PetscInt(sz);
  for (int i=0; i < sz; ++i)
    r[i] = PetscInt(rows[i]);
  CALL(VecSetValues(b,s,r,vals,ADD_VALUES));
}

void LinAlg::add_to_matrix(int sz, long* rows, double* vals) {
  PetscInt r[sz];
  PetscInt s = PetscInt(sz);
  for (int i=0; i < sz; ++i)
    r[i] = PetscInt(rows[i]);
  CALL(MatSetValues(A,s,r,s,r,vals,ADD_VALUES));
}

void LinAlg::zero_vector_rows(int sz, long* rows) {
  PetscInt r[sz];
  PetscInt s = PetscInt(sz);
  for (int i=0; i < sz; ++i)
    r[i] = PetscInt(rows[i]);
  PetscScalar vals[sz];
  for (int i=0; i < sz; ++i)
    vals[i] = 0.0;
  CALL(VecSetValues(b,s,r,vals,INSERT_VALUES));
}

void LinAlg::diag_matrix_rows(int sz, long* rows) {
  PetscInt r[sz];
  PetscInt s = PetscInt(sz);
  for (int i=0; i < sz; ++i)
    r[i] = PetscInt(rows[i]);
  CALL(MatZeroRows(A,s,r,1.0,PETSC_NULL,PETSC_NULL));
}

void LinAlg::synchronize() {
  CALL(VecAssemblyBegin(b));
  CALL(VecAssemblyEnd(b));
  CALL(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
  CALL(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
}

void LinAlg::solve() {
  double t0 = time();
  CALL(KSPSetOperators(solver,A,A));
  CALL(KSPSetFromOptions(solver));
  CALL(KSPSolve(solver,b,x));
  double t1 = time();
  long its;
  CALL(KSPGetIterationNumber(solver,&its));
  print("linear system solved in %f seconds",t1-t0);
  print("linear system converged in %ld iters",its);
}

void LinAlg::attach(Disc* d, const char* name) {
  PetscInt n;
  CALL(VecGetLocalSize(x,&n));
  PetscScalar* X;
  CALL(VecGetArray(x,&X));
  apf::GlobalNumbering* gn = d->get_numbering();
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodes(gn,nodes);
  apf::Mesh* m = d->get_apf_mesh();
  apf::Field* f = apf::createFieldOn(m, name, apf::SCALAR);
  int idx=0;
  for (unsigned i=0; i < nodes.getSize(); ++i) {
    int node = nodes[i].node;
    apf::MeshEntity* ent = nodes[i].entity;
    if (m->isOwned(nodes[i].entity))
      apf::setScalar(f,ent,node,X[idx++]);
  }
  apf::synchronize(f);
  CALL(VecRestoreArray(x,&X));
}

}
