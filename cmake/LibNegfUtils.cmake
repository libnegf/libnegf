

# Loads global build settings (either from config.cmake or from user defined file)
#
macro (libnegf_load_build_settings)

  if(NOT DEFINED BUILD_CONFIG_FILE)
	  if(DEFINED ENV{LIBNEGF_BUILD_CONFIG_FILE}
			  AND NOT "$ENV{LIBNEGF_BUILD_CONFIG_FILE}" STREQUAL "")
		  set(BUILD_CONFIG_FILE "$ENV{LIBNEGF_BUILD_CONFIG_FILE}")
    else()
      set(BUILD_CONFIG_FILE "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake")
    endif()
  endif()
  message(STATUS "Reading global build config file: ${BUILD_CONFIG_FILE}")
  include(${BUILD_CONFIG_FILE})

endmacro()
