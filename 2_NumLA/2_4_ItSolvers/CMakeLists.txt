cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)

set(SCLIB ../../SClinalg/)
set(EIGENLIB ../../eigen/)


################################
# set the project name
project(CGSolver)

# add the executable
add_executable(CGSolver CG.cpp)

target_include_directories(CGSolver PUBLIC ${SCLIB} ${EIGENLIB})
################################
################################
# set the project name
project(GMRESSolver)

# add the executable
add_executable(GMRESSolver GMRES.cpp)

target_include_directories(GMRESSolver PUBLIC ${SCLIB} ${EIGENLIB})
################################

################################
# set the project name
project(2Dconvdiff)

# add the executable
add_executable(2Dconvdiff 2D_convection_diffusion.cpp)

target_include_directories(2Dconvdiff PUBLIC ${SCLIB} ${EIGENLIB})
################################

