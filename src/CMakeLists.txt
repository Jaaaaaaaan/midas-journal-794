cmake_minimum_required (VERSION 3.8)
project (VascuSynth)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -O2 -ffast-math -fopenmp")

add_subdirectory(classes)
add_subdirectory(utils)

add_executable(VascuSynth VascuSynth.cpp)

target_include_directories (VascuSynth
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(VascuSynth
    PRIVATE Eigen3::Eigen classes utils
)
