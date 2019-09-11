cmake_minimum_required(VERSION 3.5)

project(libNEGF VERSION 0.3 LANGUAGES Fortran C)

set(FYPP "fypp" CACHE STRING "Fypp preprocessor")

option(WITH_MPI "Whether MPI-parallelised library should be built" FALSE)

# Note, this option does not work yet, so leave it to false
option(WITH_INELASTIC "Whether to build with inelastic scattering" FALSE)

set(MPIFX_INCLUDE_DIR "/opt/mpifx/include" CACHE PATH
  "Include directory of the MPIFX library (only if built with MPI)")

set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH
  "Installation directory for libraries")

set(INSTALL_MOD_DIR "${CMAKE_INSTALL_PREFIX}/include/libnegf/modfiles" CACHE PATH
  "Installation directory for Fortran module files")

option(BUILD_SHARED_LIBS "Whether the library should be shared" FALSE)

option(INSTALL_INCLUDE_FILES "Whether module files should be installed" TRUE)

if(WITH_MPI)
  find_package(MPI REQUIRED)
endif()

add_subdirectory(ext_sparskit)
add_subdirectory(ext_system)
add_subdirectory(src)
