set(CPP_LIST polygon.cpp io.cpp triang.cpp mesh.cpp  BET_elimitation.cpp hashtable.cpp metrics.cpp SmallestEnclosingCircle.cpp)
set(HPP_LIST ${PROJECT_SOURCE_DIR}/include/polygon.h ${PROJECT_SOURCE_DIR}/include/io.h ${PROJECT_SOURCE_DIR}/include/triang.h ${PROJECT_SOURCE_DIR}/include/mesh.h  ${PROJECT_SOURCE_DIR}/include/BET_elimitation.h ${PROJECT_SOURCE_DIR}/include/hashtable.h ${PROJECT_SOURCE_DIR}/include/metrics.h ${PROJECT_SOURCE_DIR}/include/SmallestEnclosingCircle.hpp)


add_library(polylla_files ${HPP_LIST} ${CPP_LIST} )

target_include_directories(polylla_files PUBLIC ../include)

#target_link_libraries(polylla_files detri2) 
#target_include_directories(polylla_files PRIVATE ../external/detri2)