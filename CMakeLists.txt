cmake_minimum_required(VERSION 3.18)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
include(LibNegfUtils)
libnegf_load_build_settings()

# Some setting output
message(STATUS "libNEGF WITH_MPI: " ${WITH_MPI})
message(STATUS "libNEGF WITH_GPU: " ${WITH_TRANSPORT_GPU})
#message(STATUS "libNEGF WITH_HILBERT: " ${WITH_HILBERT})

if(WITH_TRANSPORT_GPU)
 project(libNEGF VERSION 1.2.1 LANGUAGES Fortran C CUDA CXX)
else()
 project(libNEGF VERSION 1.2.1 LANGUAGES Fortran C CXX)
endif()


# disable compiler-specific extensions
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Note: enforcing a standard in Fortran compilers has a different effect that in C compilers, c.f. [cmake#22235](https://gitlab.kitware.com/cmake/cmake/-/issues/22235).


include(externalMpifx)
include(CMakePackageConfigHelpers)


# find a Python3 interpreter for FYPP
include(FindPython3)
find_package(Python3 REQUIRED COMPONENTS Interpreter)


set(CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/ext_fypp/;${CMAKE_PREFIX_PATH}")
find_program(FYPP "fypp")

include(GNUInstallDirs)


if(LAPACK_LIBRARIES AND LAPACK_LIBRARY_DIRS)
  message(STATUS "libNEGF LAPACK DIR: " ${LAPACK_LIBRARY_DIRS})
  message(STATUS "libNEGF BLAS: " ${BLAS_LIBRARY})
  message(STATUS "libNEGF LAPACK: " ${LAPACK_LIBRARY})
else()
  message(STATUS "libNEGF MKLROOT: " $ENV{MKLROOT})
  find_package(CustomBlas REQUIRED)
  find_package(CustomLapack REQUIRED)
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
set(CMAKE_INSTALL_PREFIX "/usr/local/libnegf${libNEGF_VERSION_MAJOR}.${libNEGF_VERSION_MINOR}" CACHE PATH
	"default installation path")

message(STATUS "Default installation path: " ${CMAKE_INSTALL_PREFIX})

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
