# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(Ex4)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)

# Add include directories (if you have any custom headers)
include_directories(include)

# Add your source files
set(SOURCES_A1
    SCTridiagLUSolver.h  
    Aufgabe1.cpp 
)

# Create an executable from your sources
add_executable(Aufgabe1 ${SOURCES_A1})

# # Add your source files
# set(SOURCES_A2
#     SCGaussSeidelSolver.h
#     TridiagGSSolverTest.cpp    
# )

# # Create an executable from your sources
# add_executable(TridiagGSSolverTest ${SOURCES_A2})

# # Add your source files
# set(SOURCES_A3_T1
#     SCGaussSeidelSolver.h
#     ApplicationExample_Task1.cpp    
# )

# # Create an executable from your sources
# add_executable(ApplicationExample ${SOURCES_A3_T1})

# # Add your source files
# set(SOURCES_A3_T2
#     SCGaussSeidelSolver.h
#     SCTridiagLUSolver.h
#     ApplicationExample_Task2.cpp    
# )

# # Create an executable from your sources
# add_executable(ApplicationExampleTiming ${SOURCES_A3_T2})