# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(Ex2)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)

# Add include directories (if you have any custom headers)
include_directories(include)

# Add your source files
set(SOURCES
    A2.cpp    
    HermiteBandMatrix.h
)

# Create an executable from your sources
add_executable(A2 ${SOURCES})

# Add your source files
set(SOURCES
    A3.cpp    
    HermiteBandMatrix.h
)

# Create an executable from your sources
add_executable(A3 ${SOURCES})
