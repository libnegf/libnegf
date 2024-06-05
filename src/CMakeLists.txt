set(sources-fpp
  clock.F90
  ln_extract.F90
  complexbands.F90
  input_output.F90
  ln_cache.F90
  lib_param.F90
  ln_enums.F90
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
  iterative.F90
  ln_constants.F90
  mat_def.F90
  population.F90
  energy_mesh.F90
  equiv_kpoints.F90
  interactions.F90
  ln_elastic.F90
  ln_inelastic.F90
  elphdd.F90
  elphdb.F90
  elphds.F90
  elphinel.F90
  self_energy.F90
  scba.F90
  sparskit/skit_blassm.F90
  sparskit/skit_formats.F90
  sparskit/skit_module.F90
  sparskit/skit_unary.F90)

if(WITH_HILBERT)
  list(APPEND sources-fpp transform.F90)
endif()

if(WITH_TRANSPORT_GPU)
  list(APPEND sources-fpp cudautils.F90)
  list(APPEND sources-fpp iterative_gpu.F90)
  set(sources-cuda gpuroutines.cu)
else()
  list(APPEND sources-fpp iterative_cpu.F90)
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


if(WITH_TRANSPORT_GPU)
  list(APPEND fyppdefs -DGPU)
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


  if(WITH_TRANSPORT_GPU)
     add_library(negf ${sources-fpp-f90} ${sources-f90} ${sources-cuda})
  else()
     add_library(negf ${sources-fpp-f90} ${sources-f90})
  endif()

target_sources(negf PRIVATE $<TARGET_OBJECTS:syscalls_objlib>)
target_include_directories(negf PUBLIC
  $<BUILD_INTERFACE:$<TARGET_PROPERTY:syscalls_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)

set(BUILD_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(negf PROPERTIES Fortran_MODULE_DIRECTORY ${BUILD_MOD_DIR})
set_property(TARGET negf PROPERTY CUDA_ARCHITECTURES OFF)

if(WITH_OMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(negf PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(WITH_MPI)
  find_package(MPI REQUIRED)
  # MPI::MPI_Fortran cannot be made PUBLIC because the compile options of this
  # target may cause errors or warnings when passed to a C compiler.
  # Consequently, we need to make these compile options conditional on the
  # language used by a dependee.
  target_link_libraries(negf PRIVATE MpiFx::MpiFx MPI::MPI_Fortran)
  target_include_directories(
    negf INTERFACE $<$<COMPILE_LANGUAGE:C,CXX>:${MPI_C_INCLUDE_DIRS}>
  )
  target_compile_definitions(
    negf INTERFACE $<$<COMPILE_LANGUAGE:C,CXX>:${MPI_C_COMPILE_DEFINITIONS}>
  )
  target_compile_options(
    negf INTERFACE $<$<COMPILE_LANGUAGE:C,CXX>:${MPI_C_COMPILE_OPTIONS}>
  )
  target_include_directories(
    negf INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_INCLUDE_DIRS}>)
  target_compile_definitions(
    negf
    INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_COMPILE_DEFINITIONS}>
  )
  target_compile_options(
    negf INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:${MPI_Fortran_COMPILE_OPTIONS}>
  )
  target_link_options(negf INTERFACE ${MPI_C_LINK_FLAGS})
  target_link_libraries(negf INTERFACE ${MPI_C_LIBRARIES})
endif()

if(WITH_TRANSPORT_GPU)
  find_package(CUDAToolkit REQUIRED)
  target_link_libraries(negf PRIVATE CUDA::cudart CUDA::cusolver CUDA::cublas)
endif()

target_link_libraries(negf PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

target_include_directories(negf PUBLIC
  $<BUILD_INTERFACE:${BUILD_MOD_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MOD_DIR}>)

install(TARGETS negf
  EXPORT negf-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${BUILD_MOD_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MOD_DIR})
endif()
