set(HEADERS
control.hpp
input.hpp
mesh.hpp
linear_algebra.hpp
tau.hpp
integrator.hpp
solver.hpp
)

set(SOURCES
control.cpp
mesh.cpp
linear_algebra.cpp
tau.cpp
integrator.cpp
solver.cpp
)

add_library(vms ${SOURCES})
target_include_directories(vms PUBLIC ${PETSC_INCLUDE_DIR})
target_link_libraries(vms PUBLIC ${PETSC_LIBRARIES})
target_link_libraries(vms PUBLIC SCOREC::core)

add_executable(1D_global_qoi 1D_global_qoi.cpp)
target_link_libraries(1D_global_qoi vms)
