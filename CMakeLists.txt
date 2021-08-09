cmake_minimum_required(VERSION 3.17)
project(BiLAD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(HEADER_LIST "${BiLAD_SOURCE_DIR}/include/bilad.h")
add_library(BiLAD src/utils.cpp src/biweight_dijkstra.cpp src/bilad.cpp src/exact_bilad.cpp src/PSQSR.cpp src/new_bilad.cpp ${HEADER_LIST})
target_include_directories(BiLAD PUBLIC ${BiLAD_SOURCE_DIR}/include)

# test
option(BUILD_TEST "build the test" ON)
if(BUILD_TEST)
    message(STATUS "setting up googletest and building tests")
    enable_testing()
    add_subdirectory(test)
endif()