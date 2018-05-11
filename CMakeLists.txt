cmake_minimum_required(VERSION 2.6)

Project(fortran_project Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

add_subdirectory(random_module)

add_executable(lattice_run src/config_m.f90 src/lattice_model.f90 src/lattice_run.f90)
target_link_libraries(lattice_run threefry_module)
