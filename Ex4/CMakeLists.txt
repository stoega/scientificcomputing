# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(Ex4)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(../eigen/)

# Create an executable from your sources
add_executable(A1 Aufgabe1.cpp)

# set the project name
project(Aufgabe2)

add_executable(A2 Aufgabe2.cpp)

# Find LAPACK
find_package(LAPACK REQUIRED)

# Link LAPACK to the executable
target_link_libraries(A2 ${LAPACK_LIBRARIES})
# target_include_directories(A2 PUBLIC ${EIGENLIB})

# set the project name
project(Aufgabe2Eigen)
include_directories(${CMAKE_SOURCE_DIR})
include_directories(../eigen/)
# add the executable
add_executable(A2_eigen Aufgabe2_eigen.cpp)
target_include_directories(A2_eigen PUBLIC ${EIGENLIB})

# Add the executable
add_executable(LAPACKExample lapack_test.cpp)
target_link_libraries(LAPACKExample ${LAPACK_LIBRARIES})