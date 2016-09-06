set(HEADERS
control.hpp
)

set(SOURCES
control.cpp
)

add_library(vms ${SOURCES})
target_include_directories(vms PUBLIC ${PETSC_INCLUDE_DIR})
target_link_libraries(vms PUBLIC ${PETSC_LIBRARIES})
target_link_libraries(vms PUBLIC SCOREC::core)

add_executable(hello hello.cpp)
target_link_libraries(hello vms)
