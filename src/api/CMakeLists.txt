set(negf-interface-headers)
set(all-header-deps libnegf_api.f90 bind_fortran)

add_custom_target(negf-c-binding ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/libnegf.h)
set(c-header-binding-deps binding/map_negf_c.txt binding/begin_c.txt binding/end_c.txt)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libnegf.h
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${all-header-deps} ${CMAKE_CURRENT_SOURCE_DIR}/${c-header-binding-deps}
    COMMENT "Generating C header"
    COMMAND ./bind_fortran -f c -m binding/map_negf_c.txt -b binding/begin_c.txt -n -t -e binding/end_c.txt libnegf_api.f90 > ${CMAKE_CURRENT_BINARY_DIR}/libnegf.h
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
list(APPEND negf-interface-headers ${CMAKE_CURRENT_BINARY_DIR}/libnegf.h)

add_custom_target(negf-cpp-binding ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/libnegf.hpp)
set(cpp-header-binding-deps binding/map_negf_cpp_wrapped.txt binding/begin_cpp_wrapped.txt binding/end_cpp_wrapped.txt)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libnegf.hpp
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${all-header-deps} ${CMAKE_CURRENT_SOURCE_DIR}/${cpp-header-binding-deps}
    COMMENT "Generating C++ header"
    COMMAND ./bind_fortran -f c++ -m binding/map_negf_cpp_wrapped.txt -b binding/begin_cpp_wrapped.txt -t -e binding/end_cpp_wrapped.txt libnegf_api.f90 > ${CMAKE_CURRENT_BINARY_DIR}/libnegf.hpp
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
list(APPEND negf-interface-headers ${CMAKE_CURRENT_BINARY_DIR}/libnegf.hpp)

add_custom_target(negf-all-headers ALL DEPENDS negf-c-binding negf-cpp-binding ${CMAKE_CURRENT_BINARY_DIR}/lnParams.h)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lnParams.h
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/lnParams.h
    COMMENT "Copying lnParams.h"
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/lnParams.h ${CMAKE_CURRENT_BINARY_DIR}/lnParams.h
)
list(APPEND negf-interface-headers ${CMAKE_CURRENT_BINARY_DIR}/lnParams.h)

if(INSTALL_INCLUDE_FILES)
    install(FILES ${negf-interface-headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()
