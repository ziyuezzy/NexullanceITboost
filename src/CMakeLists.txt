cmake_minimum_required(VERSION 3.10)
project(Nexullance_IT_cpp)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler flags
# For debugging
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -fexceptions")
# For release
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math -Wall -pedantic -fPIC")

set(Python3_ROOT_DIR /usr/bin/python3.12)
# set(Python3_EXECUTABLE /usr/bin/python3.12)
# set(Python3_LIBRARIES /usr/bin/python3.12.so)
set(Python3_INCLUDE_DIRS /usr/include/python3.12)
find_package(Python3 3.12 REQUIRED COMPONENTS Interpreter Development) 
set(pybind11_DIR /users/ziyzhang/.local/lib/python3.12/site-packages/pybind11/share/cmake/pybind11) # TODO: set up the pybind11 path
find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

# # Include directories
# include_directories(/users/ziyzhang/topology-research/nexullance/IT_boost/src/boost)
# # Link directories and libraries
# link_directories(/users/ziyzhang/topology-research/nexullance/IT_boost/src/boost/stage/lib)

# Find all source files in the current directory
file(GLOB SRCS *.cpp)

# Exclude main.cpp
list(REMOVE_ITEM SRCS ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

pybind11_add_module(${PROJECT_NAME} SHARED ${SRCS})


target_compile_definitions(${PROJECT_NAME} PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
# target_include_directories(${PROJECT_NAME} PRIVATE ${PYBIND11_INCLUDE_DIRS} /users/ziyzhang/topology-research/nexullance/IT_boost/src/boost)
target_link_libraries(${PROJECT_NAME} PRIVATE boost_graph Eigen3::Eigen)