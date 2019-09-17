set(sources
  sys_calls.c
  system_calls.f90)

add_library(syscalls_objlib OBJECT ${sources})

set(BUILD_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(syscalls_objlib PROPERTIES Fortran_MODULE_DIRECTORY ${BUILD_MOD_DIR})
set_target_properties(syscalls_objlib PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})

target_include_directories(syscalls_objlib PUBLIC
  $<BUILD_INTERFACE:${BUILD_MOD_DIR}>
  $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${BUILD_MOD_DIR}/ DESTINATION ${INSTALL_MOD_DIR})
endif()
