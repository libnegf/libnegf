cmake_minimum_required(VERSION 3.5)

# Allow defining cache variables as normal variables in a host project.
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    cmake_policy(SET CMP0077 NEW)
endif()

project(libNEGF VERSION 0.3 LANGUAGES Fortran C)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(externalMpifx)

set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/ext_fypp/;${CMAKE_PREFIX_PATH}")
find_program(FYPP "fypp")

set(FYPP "fypp" CACHE STRING "Fypp preprocessor")

option(WITH_MPI "Whether MPI-parallelised library should be built" FALSE)
message(STATUS "libNEGF WITH_MPI: " ${WITH_MPI})

# Note, this option does not work yet, so leave it to false
option(WITH_INELASTIC "Whether to build with inelastic scattering" FALSE)
message(STATUS "libNEGF WITH_INELASTIC: " ${WITH_INELASTIC})

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

option(LINK_LAPACK "Link LAPACK dependencies" TRUE)

if(LAPACK_LIBRARIES AND LAPACK_LIBRARY_DIRS)
  message(STATUS "libNEGF BLAS: " ${BLAS_LIBRARIES})
  message(STATUS "libNEGF LAPACK: " ${LAPACK_LIBRARIES})
  message(STATUS "libNEGF LAPACK DIR: " ${LAPACK_LIBRARY_DIRS})
else()	
  message(STATUS "libNEGF MKLROOT: " $ENV{MKLROOT})
  find_package(LAPACK REQUIRED)
  message(STATUS "libNEGF BLAS found: " ${BLAS_LIBRARIES})
  message(STATUS "libNEGF LAPACK found: " ${LAPACK_LIBRARIES})
  message(STATUS "libNEGF LAPACK DIR: " ${LAPACK_LIBRARY_DIRS})
endif()

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
