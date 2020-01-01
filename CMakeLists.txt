cmake_minimum_required(VERSION 3.5)

project(libNEGF VERSION 0.3 LANGUAGES Fortran C)

find_program(FYPP "fypp" PATHS "${CMAKE_SOURCE_DIR}/ext_fypp/bin")

set(FYPP "fypp" CACHE STRING "Fypp preprocessor")

option(WITH_MPI "Whether MPI-parallelised library should be built" FALSE)

# Note, this option does not work yet, so leave it to false
option(WITH_INELASTIC "Whether to build with inelastic scattering" FALSE)

set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH
  "Installation directory for libraries")

set(INSTALL_MOD_DIR "${CMAKE_INSTALL_PREFIX}/include/libnegf/modfiles" CACHE PATH
  "Installation directory for Fortran module files")

set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include/libnegf" CACHE PATH
  "Installation directory for C and C++ header files")

option(BUILD_SHARED_LIBS "Whether the library should be shared" FALSE)

option(INSTALL_INCLUDE_FILES "Whether module files and headers should be installed" TRUE)

option(BUILD_TESTING "Whether the tests should be built" FALSE)

if(WITH_MPI)
  find_package(MPI REQUIRED)
  set(MPIFX_INCLUDE_DIR "/usr/local/include/dftd3" CACHE PATH
    "Include directory of the MPIFX library (only if built with MPI)")
endif()

add_subdirectory(ext_sparskit)
add_subdirectory(ext_system)
add_subdirectory(src)
add_subdirectory(src/api)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()
