set(HEADERS
input.hpp
control.hpp
disc.hpp
linalg.hpp
tau.hpp
integrator.hpp
solver.hpp
estimator.hpp
size.hpp
adapater.hpp
)

set(SOURCES
control.cpp
disc.cpp
linalg.cpp
tau.cpp
integrator.cpp
solver.cpp
estimator.cpp
size.cpp
adapter.cpp
)

add_library(vms ${SOURCES})
target_include_directories(vms PUBLIC ${PETSC_INCLUDE_DIR})
target_link_libraries(vms PUBLIC ${PETSC_LIBRARIES})
target_link_libraries(vms PUBLIC SCOREC::core)

add_executable(lshape_globalqoi lshape_globalqoi.cpp)
target_link_libraries(lshape_globalqoi vms)

add_executable(lshape_squareqoi lshape_squareqoi.cpp)
target_link_libraries(lshape_squareqoi vms)
