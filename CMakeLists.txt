cmake_minimum_required(VERSION 3.0.0)
project(DelaunayPolyGenerator VERSION 0.1.0)
set (CMAKE_CXX_STANDARD 14)

add_subdirectory(detri2)
add_executable(DelaunayPolyGenerator main.cpp polygon.cpp io.cpp triang.cpp mesh.cpp delaunay.cpp BET_elimitation.cpp hashtable.cpp metrics.cpp SmallestEnclosingCircle.cpp)

target_link_libraries(DelaunayPolyGenerator detri2)



