include(ExternalProject)
find_package(Git REQUIRED)

function(find_or_build_mpifx)

  if (MPIFX_LIB AND MPIFX_INCLUDE_DIR)
    return()
  endif()

  find_library(
      MPIFX_LIB
      NAMES mpifx
      PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64
      DOC "Find system mpifx")
      # Figure out the include path and verify that at least the main module is
      # there.
  get_filename_component(MPIFX_DIRECTORY "${MPIFX_LIBRARY}" DIRECTORY)
  set(MPIFX_INCLUDE_DIR "${MPIFX_DIRECTORY}/../include")

  if (NOT MPIFX_LIB OR NOT EXISTS "${MPIFX_INCLUDE_DIR}/libmpifx_module.mod" OR FORCE_MPIFX_INSTALL)

    message("mpifx not found: adding a target to install from github")
    set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external/mpifx)
    set(MPIFX_INCLUDE_DIR ${EXTERNAL_INSTALL_LOCATION}/include PARENT_SCOPE)
    set(MPIFX_LIB "${EXTERNAL_INSTALL_LOCATION}/lib/libmpifx.a" PARENT_SCOPE)
    ExternalProject_Add(ext_mpifx
        GIT_REPOSITORY https://github.com/dftbplus/mpifx/
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
        -DCMAKE_POSITION_INDEPENDENT_CODE=${BUILD_SHARED_LIBS}
	  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
	  -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
    -DCMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}
    -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
    )

    add_library(mpifx STATIC IMPORTED)
    set_target_properties(mpifx PROPERTIES IMPORTED_LOCATION ${MPIFX_LIB})

  endif()

endfunction()
