set(sources-fpp
  clock.F90
  ln_extract.F90
  complexbands.F90
  input_output.F90
  lib_param.F90
  elph.F90
  phph.F90
  ln_precision.F90
  mpi_globals.F90
  contselfenergy.F90
  inversions.F90
  libnegf.F90
  ln_structure.F90
  outmatrix.F90
  rcm_module.F90
  distributions.F90
  ln_allocation.F90
  load.F90
  sparsekit_drv.F90
  integrations.F90
  globals.F90
  iterative_dns.F90
  ln_constants.F90
  mat_def.F90
  population.F90
  energy_mesh.F90
  interactions.F90
  elphdd.F90
  elphdb.F90
  elphds.F90)

if(WITH_INELASTIC)
  list(APPEND sources-fpp transform.F90)
endif()

set(sources-f90
  api/libnegfAPICommon.f90
  api/libnegf_api.f90)

#execute_process(COMMAND git describe OUTPUT_VARIABLE gitrevision)
set(gitrevision "000")
string(TIMESTAMP compdate "%Y-%m-%d")

set(fyppdefs -D_GITREVISION=\"${gitrevision}\" -D_COMPDATE=\"${compdate}\")
if(WITH_MPI)
  list(APPEND fyppdefs -DMPI)
endif()

set(sources-fpp-f90)
foreach(fppsrc IN LISTS sources-fpp)
  string(REGEX REPLACE "\\.F90" ".f90" f90src ${fppsrc})
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${f90src}
    COMMAND ${FYPP} ${fyppdefs} ${CMAKE_CURRENT_SOURCE_DIR}/${fppsrc}
      ${CMAKE_CURRENT_BINARY_DIR}/${f90src}
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${fppsrc}
    VERBATIM)
  list(APPEND sources-fpp-f90 ${CMAKE_CURRENT_BINARY_DIR}/${f90src})
endforeach()

add_library(negf ${sources-fpp-f90} ${sources-f90} )

target_sources(negf PRIVATE $<TARGET_OBJECTS:syscalls_objlib> $<TARGET_OBJECTS:sparskit>)
target_include_directories(negf PUBLIC
  $<BUILD_INTERFACE:$<TARGET_PROPERTY:syscalls_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)

target_link_libraries(negf PRIVATE sparskit)

set(BUILD_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(negf PROPERTIES Fortran_MODULE_DIRECTORY ${BUILD_MOD_DIR})

if(WITH_MPI)
  target_include_directories(negf PRIVATE ${MPIFX_INCLUDE_DIR})
  target_link_libraries(negf PRIVATE mpifx)
endif()

target_include_directories(negf PUBLIC
  $<BUILD_INTERFACE:${BUILD_MOD_DIR}>
  $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)

install(TARGETS negf
  EXPORT ${INSTALL_EXPORT_NAME}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR})

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${BUILD_MOD_DIR}/ DESTINATION ${INSTALL_MOD_DIR})
endif()
