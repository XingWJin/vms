#ifndef vms_linalg_hpp
#define vms_linalg_hpp

#include <petsc.h>

namespace vms {

class Disc;

class LinAlg {
  public:
    LinAlg(int n, long N);
    ~LinAlg();
    void add_to_vector(int sz, long* rows, double* vals);
    void add_to_matrix(int sz, long* rows, double* vals);
    void zero_vector_rows(int sz, long* rows);
    void diag_matrix_rows(int sz, long* rows);
    void synchronize();
    void solve();
    void attach(Disc* d, const char* n);
  private:
    Mat A;
    Vec x;
    Vec b;
    KSP solver;
};

}

#endif
