cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)
set(SCLIB ../../SClinalg/)

## avoid linker error when using XCode15, see https://developer.apple.com/forums/thread/737707
# if (APPLE)
# set(CMAKE_CXX_FLAGS "-no-pie")
# endif (APPLE)


################################
# set the project name
project(FloatingPoint)

# add the executable
add_executable(FloatingPoint floating_point_rep.cpp)

################################



################################
# set the project name
project(AlgorithmStability)

# add the executable
add_executable(AlgorithmStability gauss_stability_example.cpp)

target_include_directories(AlgorithmStability PUBLIC ${SCLIB})
################################

