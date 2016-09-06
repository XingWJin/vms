#ifndef vms_linear_algebra_hpp
#define vms_linear_algebra_hpp

#include <petsc.h>

namespace vms {

class Mesh;

class LinearAlgebra {
  public:
    LinearAlgebra(int n, long N);
    ~LinearAlgebra();
    void add_to_vector(int sz, long* rows, double* vals);
    void add_to_matrix(int sz, long* rows, double* vals);
    void zero_vector_rows(int sz, long* rows);
    void diag_matrix_rows(int sz, long* rows);
    void synchronize();
    void solve();
    void attach(Mesh* m);
  private:
    Mat A;
    Vec x;
    Vec b;
    KSP solver;
};

}

#endif
