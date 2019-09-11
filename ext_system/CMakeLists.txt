set(sources
  sys_calls.c
  system_calls.f90)

add_library(syscalls ${sources})

set(BUILD_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(syscalls PROPERTIES Fortran_MODULE_DIRECTORY ${BUILD_MOD_DIR})

target_include_directories(syscalls PUBLIC
  $<BUILD_INTERFACE:${BUILD_MOD_DIR}>
  $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)

install(TARGETS syscalls
  EXPORT ${INSTALL_EXPORT_NAME}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR})

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${BUILD_MOD_DIR}/ DESTINATION ${INSTALL_MOD_DIR})
endif()
