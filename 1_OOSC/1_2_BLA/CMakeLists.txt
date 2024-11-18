cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)

## avoid linker error when using XCode15, see https://developer.apple.com/forums/thread/737707
# if (APPLE)
# set(CMAKE_CXX_FLAGS "-no-pie")
# endif (APPLE)

set(SCLIB ../../SClinalg/)
set(EIGENLIB ../../eigen/)


################################
# set the project name
project(Vector)

# add the executable
add_executable(Vector vector.cpp)

# include directory for SCvector.h is needed to compile vector.cpp
target_include_directories(Vector PUBLIC ${SCLIB})
################################


################################
# set the project name
project(Vectorefficiency)

# add the executable
add_executable(Vectorefficiency vector_operator.cpp)

# include directory for SCvector.h and Eigen lib are needed
target_include_directories(Vectorefficiency PUBLIC ${SCLIB})
target_include_directories(Vectorefficiency PUBLIC ${EIGENLIB})
################################

################################
project(Matrix)

add_executable(Matrix matrix.cpp)

target_include_directories(Matrix PUBLIC ${SCLIB})
################################


################################
project(Matrixcasts)

add_executable(Matrixcasts matrix_cast.cpp)

target_include_directories(Matrixcasts PUBLIC ${SCLIB})
################################


################################
project(Matrixefficiency)

add_executable(Matrixefficiency matrix_operator.cpp)

target_include_directories(Matrixefficiency PUBLIC ${SCLIB})
target_include_directories(Matrixefficiency PUBLIC ${EIGENLIB})
################################
