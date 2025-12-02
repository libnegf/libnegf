# On some HPC systems (e.g., LUMI), front-end nodes do not feature GPUs or the front-end nodes may
# have an instruction set architecture different from compute nodes (e.g., on Fugaku with its Intel
# front-end and Arm64FX compute nodes). On these systems, GPU tests or all tests, respectively, have
# to be submitted with mpiexec/mpirun.
option(LIBNEGF_TEST_WITH_MPIEXEC "Always launch tests with mpiexec/mpirun" OFF)

# Test input files may be large and are therefore managed by git-lfs.
find_program(GIT_LFS git-lfs)

if(NOT GIT_LFS)
  message(FATAL_ERROR "Tests cannot be run because git-lfs is missing. "
                      "Disable testing or install git-lfs."
  )
endif()

# List of active tets. Each test indicates a subdirectory.
set(test-directories
    f90init
    c_mpi_init
    cpp_mpi_init
    c_int
    c_int_elph_deph
    c_int_file
    f90int
    f90int_file
    f90_transmission
    f90elph_deph
    f90read_hs
    f90Si2x2
    f90Si_nin
    f90Si_nin_40pl
)

if(WITH_TRANSPORT_GPU)
  list(APPEND test-directories testCUDA)
  list(APPEND test-directories testCUDA_decimation)
endif()

# Define a function wich copies over the content of the tests. This is used to support out-of-source
# build, as the test directory contain input files which we don't want to specify.
function(transfer_test_data testname)
  if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    add_custom_target(${testname}-data-transfer)
  else()
    add_custom_target(
      ${testname}-data-transfer
      COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR} for out-of-source build."
      COMMAND cp -f ${CMAKE_CURRENT_SOURCE_DIR}/* ${CMAKE_CURRENT_BINARY_DIR}/
    )
  endif()
  add_dependencies(${testname} ${testname}-data-transfer)
endfunction()

function(unzip_test_data testname)
  find_program(XZ NAMES xz)
  if(NOT XZ)
    message(FATAL_ERROR "cannot run test '${testname}' because xz executable is missing")
  endif()

  if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    add_custom_target(${testname}-data-unzip)
  else()
    add_custom_target(
      ${testname}-data-unzip
      COMMENT "unzip ${CMAKE_CURRENT_SOURCE_DIR} for out-of-source build."
      COMMAND tar -xJf ${CMAKE_CURRENT_SOURCE_DIR}/*.tar.xz -C ${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()
  add_dependencies(${testname}-data-unzip ${testname}-data-transfer)
  add_dependencies(${testname} ${testname}-data-unzip)
endfunction()

function(setup_test testname source_main)
  set(options GPUONLY MPI)
  set(multi_value_args SOURCES LABELS)
  cmake_parse_arguments(PARSE_ARGV 2 LIBNEGF_TEST "${options}" "" "${multi_value_args}")

  if(LIBNEGF_TEST_UNPARSED_ARGUMENTS)
    message(WARNING "could not parse the following arguments: ${LIBNEGF_TEST_UNPARSED_ARGUMENTS}")
  endif()
  if(LIBNEGF_TEST_KEYWORDS_MISSING_VALUES)
    message(
      WARNING
        "the following keyword arguments are missing values: ${LIBNEGF_TEST_KEYWORDS_MISSING_VALUES}"
    )
  endif()

  list(PREPEND LIBNEGF_TEST_SOURCES ${source_main})
  if(LIBNEGF_TEST_GPUONLY)
    list(APPEND LIBNEGF_TEST_LABELS gpu-only)
  endif()
  if(LIBNEGF_TEST_MPI)
    list(APPEND LIBNEGF_TEST_LABELS mpi)
  endif()

  # detect implementation language
  if(source_main MATCHES "[.]c$")
    set(LIBNEGF_TEST_LANG C)
  elseif(source_main MATCHES "[.]cpp$")
    set(LIBNEGF_TEST_LANG CXX)
  elseif(source_main MATCHES "[.](f|F)90$")
    set(LIBNEGF_TEST_LANG Fortran)
  else()
    message(SEND_ERROR "unknown file extension of test source ${source_main}")
  endif()

  message(DEBUG "add test: ${testname} ${LIBNEGF_TEST_SOURCES} ${LIBNEGF_TEST_LABELS}")

  # set up test compilation
  add_executable(${testname} ${LIBNEGF_TEST_SOURCES})

  if(LIBNEGF_TEST_LANG STREQUAL C)
    target_link_libraries(${testname} MPI::MPI_C)
    target_include_directories(${testname} PRIVATE ${CMAKE_BINARY_DIR}/src/api)
  elseif(LIBNEGF_TEST_LANG STREQUAL CXX)
    target_link_libraries(${testname} MPI::MPI_CXX)
    target_include_directories(${testname} PRIVATE ${CMAKE_BINARY_DIR}/src/api)
  elseif(LIBNEGF_TEST_LANG STREQUAL Fortran)
    target_link_libraries(${testname} MPI::MPI_Fortran)
    target_include_directories(${testname} PRIVATE ${BUILD_MOD_DIR})
  endif()
  target_link_libraries(${testname} negf)
  target_link_libraries(${testname} LAPACK::LAPACK)

  # determine how to execute the test
  if(${LIBNEGF_TEST_MPI})
    set(num_procs ${MPIEXEC_MAX_NUMPROCS})
    set(cmd
        # When running on a user's PC, it is possible to thrash the machine by having one MPI task
        # per core and one OpenMP thread per task and per core (the default). The restriction below
        # aims to avoid this scenario.
        #
        # env is wrapping the call to mpiexec/mpirun to avoid the job name `env` when
        # ${MPIEXEC_EXECUTABLE} is a batch scheduler.
        env
        OMP_NUM_THREADS=1
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG}
        ${num_procs}
        ${MPIEXEC_PREFLAGS}
        ./${testname}
        ${MPIEXEC_POSTFLAGS}
    )
  elseif(${LIBNEGF_TEST_WITH_MPIEXEC} AND NOT ${LIBNEGF_TEST_MPI})
    set(num_procs 1)
    set(cmd ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${num_procs} ${MPIEXEC_PREFLAGS}
            ./${testname} ${MPIEXEC_POSTFLAGS}
    )
  else()
    set(cmd ${testname})
  endif()

  add_test(NAME ${testname} COMMAND ${cmd})
  set_property(TEST ${testname} PROPERTY LABELS ${LIBNEGF_TEST_LABELS})
  set_property(TEST ${testname} PROPERTY PROCESSORS ${num_procs})

  transfer_test_data(${testname})
endfunction()

foreach(test-directory IN LISTS test-directories)
  add_subdirectory(${test-directory})
endforeach()
