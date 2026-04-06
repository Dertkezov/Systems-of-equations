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

add_executable(solution main.cpp)
target_link_libraries(solution PRIVATE matrix_lib)

add_executable(matrix matrix.cpp)
target_link_libraries(matrix PRIVATE matrix_lib)

add_executable(matrix_operations matrix_operations.cpp)
target_link_libraries(matrix_operations PRIVATE matrix_lib)

add_executable(timer timer.cpp)
target_link_libraries(timer PRIVATE matrix_lib)

add_executable(solver solver.cpp)
target_link_libraries(solver PRIVATE matrix_lib)

add_executable(compare compare.cpp)
target_link_libraries(compare PRIVATE matrix_lib)

add_executable(compare_symmetric compare_symmetric.cpp)
target_link_libraries(compare_symmetric PRIVATE matrix_lib)

add_executable(symmetric_solver symmetric_solver.cpp)
target_link_libraries(symmetric_solver PRIVATE matrix_lib)

add_executable(compare_ellipse_matrix compare_ellipse_matrix.cpp)
target_link_libraries(compare_ellipse_matrix PRIVATE matrix_lib)