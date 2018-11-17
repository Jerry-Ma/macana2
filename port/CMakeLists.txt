project(port_test LANGUAGES CXX C)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(port_test EXCLUDE_FROM_ALL)
set_target_properties(port_test
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
target_sources(port_test
    PRIVATE
        port_generic_curvefit.cpp
        test/test_port_generic_curvefit.cpp
    )
target_include_directories(port_test
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
    )
target_link_libraries(port_test
    PRIVATE
        gtest gtest_main
        Eigen3::Eigen
    )
add_custom_command(TARGET port_test
    POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
    )
add_dependencies(check port_test)
gtest_discover_tests(port_test TEST_PREFIX "port::")
