cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)

set(SCLIB ../../SClinalg/)
set(EIGENLIB ../../eigen/)

################################
# set the project name
project(QR)

# add the executable
add_executable(QR QR.cpp)
target_include_directories(QR PUBLIC ${SCLIB})
target_include_directories(QR PUBLIC ${EIGENLIB})

################################



################################
# set the project name
project(QRShift)

# add the executable
add_executable(QRShift QR_shift.cpp)
target_include_directories(QRShift PUBLIC ${SCLIB})
target_include_directories(QRShift PUBLIC ${EIGENLIB})

################################


################################
# set the project name
project(WaveEquation)

# add the executable
add_executable(WaveEquation waveequation.cpp)
target_include_directories(WaveEquation PUBLIC ${SCLIB})

################################



################################
# set the project name
project(WaveEquation2)

# add the executable
add_executable(WaveEquation2 waveequation_Eigen.cpp)
target_include_directories(WaveEquation2 PUBLIC ${EIGENLIB})
################################
