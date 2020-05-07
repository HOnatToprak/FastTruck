cmake_minimum_required(VERSION 3.16)
project(FASTTRUCK_ROOT)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
add_subdirectory(thirdparty/armadillo-9.870.2)
add_subdirectory(fasttruck)
add_subdirectory(test)

message(STATUS "ONAT ${CMAKE_CXX_FLAGS}")

