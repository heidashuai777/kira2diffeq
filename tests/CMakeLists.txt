add_executable(test_sunrise test_sunrise.cpp)
target_link_libraries(test_sunrise kira2diffeq)
target_compile_definitions(test_sunrise PRIVATE TEST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures/cli_smoke")
add_test(NAME test_sunrise COMMAND test_sunrise)

add_test(
    NAME cli_smoke
    COMMAND ${CMAKE_BINARY_DIR}/kira2diffeq-cli
        -c ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/cli_smoke
        -m ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/cli_smoke/reduction.m
        -f ${CMAKE_CURRENT_SOURCE_DIR}/fixtures/cli_smoke/masters.final
        -t cli_smoke
        -v t
        -o text
        -O ${CMAKE_CURRENT_BINARY_DIR}/cli_smoke_de.txt
)

set_tests_properties(cli_smoke PROPERTIES
    WILL_FAIL false
)
