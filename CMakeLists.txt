cmake_minimum_required(VERSION 3.16.3)
project(Polylla VERSION 0.1.0)
set (CMAKE_CXX_STANDARD 14)

add_subdirectory(external)
add_subdirectory(src)

add_executable(Polylla main.cpp)

target_link_libraries(Polylla polylla_files)
