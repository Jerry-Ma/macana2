project(port_test LANGUAGES CXX C)

FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.x
)
FetchContent_GetProperties(spdlog)
if(NOT spdlog_POPULATED)
    FetchContent_Populate(spdlog)
    add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
find_package(spdlog)
print_target_properties(spdlog::spdlog)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Ceres REQUIRED COMPONENTS C++11)
find_package( OpenCV REQUIRED )

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
        spdlog::spdlog
        gtest gtest_main gmock
        ceres
        Eigen3::Eigen
        ${OpenCV_LIBS}
    )
add_custom_command(TARGET port_test
    POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
    )
add_dependencies(check port_test)
gtest_discover_tests(port_test TEST_PREFIX "port::")
