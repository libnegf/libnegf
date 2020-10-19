function(find_or_build_mpifx)

  # If MpiFx had been detected already, done.
  if(TARGET MpiFx::MpiFx)
    return()
  endif()
  
  # If user has set up library explicitely, take that choice withouth further checks.
  set(MPIFX_LIBRARY "" CACHE STRING "MpiFx library path")
  set(MPIFX_INCLUDE_DIRECTORY "" CACHE STRING "MpiFx include path")
  if (MPIFX_LIBRARY AND MPIFX_INCLUDE_DIRECTORY)
    message(STATUS "Customized MpiFx library: ${MPIFX_LIBRARY}")
    message(STATUS "Customized MpiFx include directory: ${MPIFX_INCLUDE_DIRECTORY}")
    add_library(MpiFx::MpiFx INTERFACE IMPORTED GLOBAL)
    target_link_libraries(MpiFx::MpiFx INTERFACE "${MPIFX_LIBRARY}")
    target_include_directories(MpiFx::MpiFx INTERFACE "${MPIFX_INCLUDE_DIRECTORY}")
    return()
  endif()

  if(NOT FORCE_MPIFX_DOWNLOAD)
    # Find system installed MpiFx package.
    find_package(MpiFx QUIET)
    if(MpiFx_FOUND)
      message(STATUS "Found installed MpiFx package")
      return()
    endif()
  endif()

  # Download package and build it as part of normal build process
  message(STATUS "Downloading MpiFx from GitHub")
  include(FetchContent)
  FetchContent_Declare(FetchedMpiFx
    GIT_REPOSITORY https://github.com/aradi/mpifx/
    GIT_TAG cmake)
  FetchContent_MakeAvailable(FetchedMpiFx)
  add_library(MpiFx::MpiFx INTERFACE IMPORTED GLOBAL)
  target_link_libraries(MpiFx::MpiFx INTERFACE MpiFx)

endfunction()
