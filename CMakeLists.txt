cmake_minimum_required(VERSION 2.6)

Project(fortran_project Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS_INIT "-Wall -std=c99")
elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS_INIT "-std=c99")
endif()

add_subdirectory(random_module)

add_executable(lattice_run src/config_m.f90 src/lattice_model.f90 src/lattice_run.f90)
target_link_libraries(lattice_run xoroshiro128plus_m)

add_executable(ssa_run src/config_m.f90 src/ssa_model_m.f90 src/ssa_run.f90)
