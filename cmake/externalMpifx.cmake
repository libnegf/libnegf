include(ExternalProject)
find_package(Git REQUIRED)

function(find_or_build_mpifx)

  find_library(
      MPIFX_LIBRARY
      NAMES mpifx
      PATHS /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64
      DOC "Find system mpifx")
      # Figure out the include path and verify that at least the main module is
      # there.
  get_filename_component(MPIFX_DIRECTORY "${MPIFX_LIBRARY}" DIRECTORY)
  set(MPIFX_INCLUDE_DIR "${MPIFX_DIRECTORY}/../include")

  if (NOT MPIFX_LIBRARY OR NOT EXISTS "${MPIFX_INCLUDE_DIR}/libmpifx_module.mod" OR FORCE_MPIFX_INSTALL)

    message("mpifx not found: adding a target to install from github")
    set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external/mpifx)
    ExternalProject_Add(ext_mpifx
        GIT_REPOSITORY https://github.com/dftbplus/mpifx/
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
    )

    set(MPIFX_INCLUDE_DIR ${EXTERNAL_INSTALL_LOCATION}/include PARENT_SCOPE)
    set(MPIFX_LIBRARY "${EXTERNAL_INSTALL_LOCATION}/lib/libmpifx.a" PARENT_SCOPE)
  endif()

  include_directories(${MPIFX_INCLUDE_DIR})

endfunction()
