#include <cstdarg>
#include <cstdio>

#include <PCU.h>
#include <petsc.h>

namespace vms {

void initialize() {
  MPI_Init(0,0);
  PCU_Comm_Init();
  PetscInitialize(0,0,0,0);
}

void finalize() {
  PetscFinalize();
  PCU_Comm_Free();
  MPI_Finalize();
}

double time()
{
  return PCU_Time();
}

void print(const char* format, ...) {
  if (PCU_Comm_Self())
    return void();
  va_list ap;
  va_start(ap, format);
  vfprintf(stdout, format, ap);
  va_end(ap);
  printf("\n");
}

void fail(char const* format, ...) {
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  printf("\n");
  abort();
}

}
