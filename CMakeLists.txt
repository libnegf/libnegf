cmake_minimum_required(VERSION 3.5)

project(libNEGF VERSION 0.3 LANGUAGES Fortran C)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(externalMpifx)

set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/ext_fypp/;${CMAKE_PREFIX_PATH}")
find_program(FYPP "fypp")

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

option(FORCE_MPIFX_INSTALL "Install mpifx from repository even if a version is found in the system" FALSE)

# Requirements.
find_package(BLAS REQUIRED)
message(STATUS "BLAS found:" ${BLAS_LIBRARIES})
message(STATUS "BLAS vendor:" ${BLAS_VENDOR_FOUND})
message(STATUS "BLAS linker flags:" ${BLAS_LINKER_FLAGS})
find_package(LAPACK REQUIRED)
message(STATUS "LAPACK found:" ${LAPACK_LIBRARIES})
message(STATUS "LAPACK linker flags:" ${LAPACK_LINKER_FLAGS})

if(WITH_MPI)
  find_package(MPI REQUIRED)
  find_or_build_mpifx()
endif()

# Subdirectories.
add_subdirectory(ext_sparskit)
add_subdirectory(ext_system)
add_subdirectory(src)
add_subdirectory(src/api)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()
