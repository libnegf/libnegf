# List of active tets. Each test indicates a subdirectory.
set(test-directories
        c_int
        c_int_elph_deph
        c_int_file
        f90int
        f90int_file
        f90sgf_cache
        f90elph_deph
        f90issue6
        f90block_partition
        )

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

# A function to setup C tests.
function(setup_c_test testname source)
    add_executable(${testname} ${source})
    transfer_test_data(${testname})
    target_link_libraries(${testname} negf)
    target_include_directories(${testname} PRIVATE ${CMAKE_BINARY_DIR}/src/api)
    add_test(${testname} ${testname})
endfunction()

# A function to setup F90 tests.
function(setup_f90_test testname source)
    add_executable(${testname} ${source})
    transfer_test_data(${testname})
    target_link_libraries(${testname} negf)
    target_include_directories(${testname} PRIVATE ${BUILD_MOD_DIR})
    add_test(${testname} ${testname})
endfunction()


foreach(test-directory IN LISTS test-directories)
    add_subdirectory(${test-directory})
endforeach()
