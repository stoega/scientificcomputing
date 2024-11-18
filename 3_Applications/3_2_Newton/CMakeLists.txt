cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)

set(SCLIB ../../SClinalg/)

################################
# set the project name
project(Newton)

# add the executable
add_executable(Newton newton.cpp)
target_include_directories(Newton PUBLIC ${SCLIB})

################################

################################
# set the project name
project(NewtonAutoDiff)

# add the executable
add_executable(NewtonAutoDiff newton_autodiff.cpp)
target_include_directories(NewtonAutoDiff PUBLIC ${SCLIB})

################################


################################
# set the project name
project(GaussNewton)

# add the executable
add_executable(GaussNewton fitting.cpp)
target_include_directories(GaussNewton PUBLIC ${SCLIB})

################################


