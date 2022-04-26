function(find_or_fetch_mpifx)

  # If MpiFx had been detected already, done.
  if(TARGET MpiFx::MpiFx)
    return()
  endif()

  if(NOT FORCE_MPIFX_DOWNLOAD)
    # Find system installed MpiFx package (provided CMAKE_PREFIX_PATH contains the package root)
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
    GIT_REPOSITORY https://github.com/dftbplus/mpifx/
    GIT_TAG release)
  FetchContent_MakeAvailable(FetchedMpiFx)
  add_library(MpiFx::MpiFx INTERFACE IMPORTED GLOBAL)
  target_link_libraries(MpiFx::MpiFx INTERFACE MpiFx)

endfunction()
