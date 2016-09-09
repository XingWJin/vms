find_path(PETSC_INCLUDE_DIR petsc.h HINTS ${Petsc_PREFIX}/include)
find_library(PETSC_LIBRARIES NAME petsc HINTS ${Petsc_PREFIX}/lib)

if(PETSC_LIBRARIES)
set(PETSC_LIBRARY_DIR ${Petsc_PREFIX}/lib)
set(PETSC_FOUND ON)
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args (PETSC
  "can't find petsc"
PETSC_INCLUDE_DIR PETSC_LIBRARIES PETSC_LIBRARY_DIR)
