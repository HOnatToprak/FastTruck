project(fasttruck)
add_library(${PROJECT_NAME}
        src/TSP.cpp
        src/Utils.cpp)
find_package(OpenMP REQUIRED)
target_include_directories(${PROJECT_NAME}
        PUBLIC ${PROJECT_SOURCE_DIR}/include
        PUBLIC ${CMAKE_SOURCE_DIR}/thirdparty/armadillo-9.870.2/include
        PUBLIC ${OpenMP_CXX_INCLUDE_DIRS}
)



