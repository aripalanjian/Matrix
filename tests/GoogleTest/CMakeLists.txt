cmake_minimum_required(VERSION 3.31)
# 'test_g' is the subproject name
project(GoogleTest)

set(CMAKE_CXX_STANDARD 20)

# 'lib' is the folder with Google Test sources
add_subdirectory(lib build)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# 'run_gtests' is the target name
# 'operator_test.cpp' are source files with tests
add_executable(run_gtests
        operator_test.cpp
        MatrixOperatorTest.cpp)
target_link_libraries(run_gtests gtest gtest_main)


