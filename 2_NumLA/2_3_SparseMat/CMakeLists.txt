cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)

set(SCLIB ../../SClinalg/)
set(EIGENLIB ../../eigen/)


################################
# set the project name
project(SparseMatrix)

# add the executable
add_executable(SparseMatrix sparsematrixtest.cpp)

target_include_directories(SparseMatrix PUBLIC ${SCLIB} ${EIGENLIB})
################################

################################
# set the project name
project(Poisson)

# add the executable
add_executable(Poisson Poisson.cpp)
target_include_directories(Poisson PUBLIC ${SCLIB})
################################


################################
# set the project name
project(SparseSolver)

# add the executable
add_executable(SparseSolver sparsesolver_Eigen.cpp)

target_include_directories(SparseSolver PUBLIC ${SCLIB} ${EIGENLIB})
################################

################################
# set the project name
project(LU_Fill_In)

# add the executable
add_executable(LU_Fill_In LU_fill_in.cpp)

target_include_directories(LU_Fill_In PUBLIC ${SCLIB} ${EIGENLIB})
################################
