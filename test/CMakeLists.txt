project (macana2_test LANGUAGES CXX C)

add_executable(macana2_test EXCLUDE_FROM_ALL)
set_target_properties(macana2_test
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
target_sources(macana2_test
    PRIVATE
        AnalParamsTest.cpp
        GaussFitTest.cpp
        MapTest.cpp
    )
target_link_libraries(macana2_test
    PRIVATE
        macana2
        gtest gtest_main gmock
    )
add_custom_command(TARGET macana2_test
    POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
    )
add_dependencies(check macana2_test)
gtest_discover_tests(macana2_test TEST_PREFIX "macana2::")
