project(test)
add_executable(${PROJECT_NAME} main.cpp)

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME}
        fasttruck
        armadillo
        OpenMP::OpenMP_CXX
)
