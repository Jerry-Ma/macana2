project(port_test LANGUAGES CXX C)

FetchContentHelper(spdlog https://github.com/gabime/spdlog.git v1.x)
FetchContentHelper(fmt https://github.com/fmtlib/fmt.git 5.2.1)

print_target_properties(spdlog::spdlog)
print_target_properties(fmt::fmt)

find_package(OpenCV)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_COMPILE_DEFINITIONS "EIGEN_NO_MALLOC")
find_package(Ceres REQUIRED COMPONENTS C++11)

add_library(port STATIC)
target_compile_definitions(port PUBLIC _USE_MATH_DEFINES)
target_sources(port
    PRIVATE
        port_generic_curvefit.cpp
        port_timestream_sensitivity.cpp
    )
target_include_directories(port
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
    )
target_link_libraries(port
    PUBLIC
        spdlog::spdlog
        fmt::fmt
        Eigen3::Eigen
        ceres
    )

add_executable(port_test EXCLUDE_FROM_ALL)
set_target_properties(port_test
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
target_sources(port_test
    PRIVATE
        test/test_port_generic_curvefit.cpp
        test/test_port_timestream_sensitivity.cpp
        test/test_main.cpp
    )
target_include_directories(port_test
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
    )
target_link_libraries(port_test
    PRIVATE
        port
        gtest gmock benchmark
        ${OpenCV_LIBS}
    )
add_custom_command(TARGET port_test
    POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
    )
add_dependencies(check port_test)
gtest_discover_tests(port_test TEST_PREFIX "port::")
