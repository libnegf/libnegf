include(ExternalProject)
find_package(Git REQUIRED)

function(find_or_build_mpifx)

  # mpifx provided as target by top project.
  if (TARGET mpifx)
    message(STATUS "MPIFX library: ${MPIFX_LIB}")
    message(STATUS "MPIFX includes: ${MPIFX_INCLUDE_DIR}")
    return()
  endif()

  # mpifx provided as path to existing library.
  set(MPIFX_LIB "" CACHE STRING "mpifx library path")
  set(MPIFX_INCLUDE_DIR "" CACHE STRING "mpifx include path")
  if (MPIFX_LIB AND MPIFX_INCLUDE_DIR)
    message(STATUS "MPIFX library: ${MPIFX_LIB}")
    message(STATUS "MPIFX includes: ${MPIFX_INCLUDE_DIR}")
    add_library(mpifx UNKNOWN IMPORTED)
    set_property(TARGET mpifx PROPERTY IMPORTED_LOCATION ${MPIFX_LIB})
    return()
  endif()

  # mpifx searched as system library.
  find_library(
      MPIFX_LIB_FIND
      NAMES mpifx
      PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64
      DOC "Find system mpifx")
  find_path(
      MPIFX_INCLUDE_DIR_FIND
      NAMES libmpifx_module.mod
      PATHS /usr/lib/mpifx /usr/lib64/mpifx /usr/local/ /usr/local/mpifx
      PATH_SUFFIXES include includes modfiles modules
  )
  set(MPIFX_LIB "${MPIFX_LIB_FIND}")
  set(MPIFX_INCLUDE_DIR "${MPIFX_INCLUDE_DIR_FIND}")
  if (MPIFX_LIB)
    message(STATUS "MPIFX found: ${MPIFX_LIB}")
  endif()
  if (MPIFX_INCLUDE_DIR)
    message(STATUS "MPIFX include found: ${MPIFX_INCLUDE_DIR}")
  endif()

  add_library(mpifx UNKNOWN IMPORTED)

  # Everything failed: pull and compile from github.
  if (NOT MPIFX_LIB OR NOT MPIFX_INCLUDE_DIR OR FORCE_MPIFX_INSTALL)
    message(STATUS "MPIFX: installing from https://github.com/dftbplus/mpifx/")
    set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external/mpifx)
    ExternalProject_Add(ext_mpifx
    GIT_REPOSITORY https://github.com/dftbplus/mpifx/
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
    -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_SHARED_LIBS}
	  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
	  -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
    -DCMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}
    -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
    )

    set(MPIFX_INCLUDE_DIR ${EXTERNAL_INSTALL_LOCATION}/include)
    set(MPIFX_LIB "${EXTERNAL_INSTALL_LOCATION}/lib/libmpifx.a")
    add_dependencies(mpifx ext_mpifx)
  endif()

  set(MPIFX_INCLUDE_DIR ${MPIFX_INCLUDE_DIR} PARENT_SCOPE)
  set_property(TARGET mpifx PROPERTY IMPORTED_LOCATION ${MPIFX_LIB})

endfunction()
