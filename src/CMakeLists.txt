set(HEADERS
input.hpp
control.hpp
disc.hpp
linalg.hpp
tau.hpp
integrator.hpp
solver.hpp
)

set(SOURCES
control.cpp
disc.cpp
linalg.cpp
tau.cpp
integrator.cpp
solver.cpp
)

add_library(vms ${SOURCES})
target_include_directories(vms PUBLIC ${PETSC_INCLUDE_DIR})
target_link_libraries(vms PUBLIC ${PETSC_LIBRARIES})
target_link_libraries(vms PUBLIC SCOREC::core)

add_executable(lshape_globalqoi_spr lshape_globalqoi_spr.cpp)
target_link_libraries(lshape_globalqoi_spr vms)
