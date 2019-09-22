# List of active tets. Each test indicates a subdirectory. 
set(test-directories 
        c_int 
        c_int_file 
        f90int 
        f90int_file
        )

# Make BLAS and LAPACK available to be linked in the tests.
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(libs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

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


foreach(test-directory IN LISTS test-directories)
    add_subdirectory(${test-directory})
endforeach()
