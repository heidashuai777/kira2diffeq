cmake_minimum_required(VERSION 3.16)
project(kira2diffeq VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Threads REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GINAC REQUIRED ginac)
pkg_check_modules(YAMLCPP REQUIRED yaml-cpp)

add_library(kira2diffeq STATIC
    src/config_reader.cpp
    src/math_parser.cpp
    src/differentiator.cpp
    src/reducer.cpp
    src/exporter.cpp
)

target_include_directories(kira2diffeq PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${GINAC_INCLUDE_DIRS}
    ${YAMLCPP_INCLUDE_DIRS}
)

target_link_libraries(kira2diffeq PUBLIC
    ${GINAC_LIBRARIES}
    ${YAMLCPP_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

set_target_properties(kira2diffeq PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(kira2diffeq-cli src/main.cpp)
target_link_libraries(kira2diffeq-cli kira2diffeq ${CMAKE_THREAD_LIBS_INIT})

enable_testing()
add_subdirectory(tests)

install(TARGETS kira2diffeq kira2diffeq-cli
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/ DESTINATION include)
