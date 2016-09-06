cmake_minimum_required(VERSION 3.0)
project(VMS LANGUAGES CXX)

set(VMS_MAJOR "1")
set(VMS_MINOR "0")
set(VMS_PATCH "0")

set(Scorec_PREFIX "" CACHE PATH "Scorec installation directory")
set(Petsc_PREFIX "" CACHE PATH "Petsc installation directory")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

find_package(Petsc REQUIRED)
find_package(SCOREC 2.0.0 REQUIRED PATHS ${Scorec_PREFIX})

add_subdirectory(src)
