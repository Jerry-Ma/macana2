
file(GLOB SRCS *.cpp)

add_executable(test ${SRCS})


target_link_libraries(test macana-core ${link_libraries} gmock_main)

add_test(NAME test COMMAND test)
