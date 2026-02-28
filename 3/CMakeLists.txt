cmake_minimum_required(VERSION 3.14)
project(MatrixBenchmark LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_TESTS "Build unit tests" ON)

add_subdirectory(src)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

add_executable(main main.cpp)
target_link_libraries(main PRIVATE matrix_lib)
