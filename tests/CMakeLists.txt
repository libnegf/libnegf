# Test input files may be large and are therefore managed by git-lfs.
find_program(GIT_LFS git-lfs)

if(NOT GIT_LFS)
  message(FATAL_ERROR
    "Tests cannot be run because git-lfs is missing. "
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
        )

if(WITH_TRANSPORT_GPU)
  list(APPEND test-directories testCUDA)
  list(APPEND test-directories testCUDA_decimation)
endif()

# Define a function wich copies over the content of the tests.
# This is used to support out-of-source build, as the test
# directory contain input files which we don't want to specify.
function(transfer_test_data testname)
    if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
        add_custom_target(${testname}-data-transfer)
    else()
        add_custom_target(
            ${testname}-data-transfer
            COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR} for out-of-source build."
            COMMAND cp -f ${CMAKE_CURRENT_SOURCE_DIR}/* ${CMAKE_CURRENT_BINARY_DIR}/)
    endif()
    add_dependencies(${testname} ${testname}-data-transfer)
endfunction()

function(unzip_test_data testname)
    if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
	  add_custom_target(${testname}-data-unzip)
        add_dependencies(${testname}-data-unzip ${testname}-data-transfer)
    else()
        add_custom_target(
            ${testname}-data-unzip
            COMMENT "unzip ${CMAKE_CURRENT_SOURCE_DIR} for out-of-source build."
		COMMAND tar -xJf ${CMAKE_CURRENT_BINARY_DIR}/*.tar.xz
		COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/*.tar)
    endif()
    add_dependencies(${testname}-data-unzip ${testname}-data-transfer)
    add_dependencies(${testname} ${testname}-data-unzip)
endfunction()

# A function to setup C++ tests.
function(setup_cpp_test testname source)
    add_executable(${testname} ${source})
    transfer_test_data(${testname})
    target_link_libraries(${testname} MPI::MPI_CXX)
    target_link_libraries(${testname} negf)
    target_link_libraries(${testname} LAPACK::LAPACK)
    target_include_directories(${testname} PRIVATE ${CMAKE_BINARY_DIR}/src/api)
    add_test(${testname} ${testname})
endfunction()

# A function to setup C tests.
function(setup_c_test testname source)
    add_executable(${testname} ${source} ${ARGN})
    transfer_test_data(${testname})
    target_link_libraries(${testname} MPI::MPI_C)
    target_link_libraries(${testname} negf)
    target_link_libraries(${testname} LAPACK::LAPACK)
    target_include_directories(${testname} PRIVATE ${CMAKE_BINARY_DIR}/src/api)
    add_test(${testname} ${testname})
endfunction()

# A function to setup F90 tests.
function(setup_f90_test testname source)
    add_executable(${testname} ${source} ${ARGN})
    transfer_test_data(${testname})
    target_link_libraries(${testname} MPI::MPI_Fortran)
    target_link_libraries(${testname} negf)
    target_link_libraries(${testname} LAPACK::LAPACK)
    target_include_directories(${testname} PRIVATE ${BUILD_MOD_DIR})
    add_test(${testname} ${testname})
endfunction()


foreach(test-directory IN LISTS test-directories)
    add_subdirectory(${test-directory})
endforeach()
