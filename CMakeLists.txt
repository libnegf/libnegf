cmake_minimum_required(VERSION 3.16)

project(libNEGF VERSION 1.0.3 LANGUAGES Fortran C)

# Optional parameters
# Set environmnent variables FC and CC; module load should do this

# MPI options
option(WITH_MPI "Whether MPI-parallelised library should be built" FALSE)
message(STATUS "libNEGF WITH_MPI: " ${WITH_MPI})

# Note, this option does not work yet, so leave it to false
option(WITH_INELASTIC "Whether to build with inelastic scattering" FALSE)
message(STATUS "libNEGF WITH_INELASTIC: " ${WITH_INELASTIC})

# shared library .so build 
option(BUILD_SHARED_LIBS "Whether the library should be shared" FALSE)

# Whether include modules should be added to installation
option(INSTALL_INCLUDE_FILES "Whether module files and headers should be installed" TRUE)

# Whether the tests should be compiled
option(BUILD_TESTING "Whether the tests should be built" FALSE)

# Whether mpifx should be downloaded 
option(FORCE_MPIFX_DOWNLOAD
  "Force mpifx download from repository (do not search for installed package) " FALSE)

set(MKL_LIBDIR "/usr/pack/intel_mkl-11.1-ma/mkl/lib/intel64" CACHE STRING
	"Directory where to look for mkl libraries"	)

set(LAPACK_LIBRARY_DIR "${MKL_LIBDIR};/usr/lib/x86_64-linux-gnu"
	CACHE STRING "Paths to scan when looking for the BLAS/LAPACK libraries")

set(LAPACK_LIBRARY "mkl_intel_lp64;mkl_gnu_thread;mkl_core;-lgomp" CACHE STRING "BLAS/LAPACK libraries")
     
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(externalMpifx)
include(CMakePackageConfigHelpers)

set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/ext_fypp/;${CMAKE_PREFIX_PATH}")
find_program(FYPP "fypp")

include(GNUInstallDirs)

set(INSTALL_INCLUDE_DIR "libnegf" CACHE PATH
  "Installation directory for C and C++ header files (within standard include folder)")

set(INSTALL_MOD_DIR "${INSTALL_INCLUDE_DIR}/modfiles" CACHE PATH
  "Installation directory for Fortran module files (within standard include folder)")


if(LAPACK_LIBRARIES AND LAPACK_LIBRARY_DIRS)
  message(STATUS "libNEGF BLAS: " ${BLAS_LIBRARIES})
  message(STATUS "libNEGF LAPACK: " ${LAPACK_LIBRARIES})
  message(STATUS "libNEGF LAPACK DIR: " ${LAPACK_LIBRARY_DIRS})
else()
  message(STATUS "libNEGF MKLROOT: " $ENV{MKLROOT})
  find_package(CustomLapack REQUIRED)
  message(STATUS "libNEGF BLAS found: " ${BLAS_LIBRARIES})
  message(STATUS "libNEGF LAPACK found: " ${LAPACK_LIBRARIES})
  message(STATUS "libNEGF LAPACK DIR: " ${LAPACK_LIBRARY_DIRS})
endif()

if(WITH_MPI)
  find_package(MPI REQUIRED)
  find_or_fetch_mpifx()
endif()

# Subdirectories.
add_subdirectory(ext_system)
add_subdirectory(src)
add_subdirectory(src/api)
add_subdirectory(src/sparskit)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()


#
# Installation
#
add_library(Negf INTERFACE)
target_link_libraries(Negf INTERFACE negf)
install(TARGETS Negf EXPORT negf-targets)

install(EXPORT negf-targets
  FILE negf-targets.cmake
  NAMESPACE Negf::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/negf"
  EXPORT_LINK_INTERFACE_LIBRARIES)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/utils/export/negf-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/negf-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/negf)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/negf-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/negf-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/negf-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/negf)
