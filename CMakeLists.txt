cmake_minimum_required(VERSION 3.31)
project(Neo VERSION 0.1)

set(CMAKE_CXX_STANDARD 20)

add_library(Neo INTERFACE)
target_link_libraries(Neo INTERFACE include)

enable_testing()
add_subdirectory(tests)
