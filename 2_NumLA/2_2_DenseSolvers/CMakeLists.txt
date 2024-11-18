cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)

set(SCLIB ../../SClinalg/)
set(EIGENLIB ../../eigen/)

option( USE_MKL "enable MKL" off)
option( USE_LAPACK "enable LAPACK without MKL" off)

if (USE_MKL)
    # preferred: edit .vscode/settings.json
    # quick and dirty:
    # manually setting MKL_INCLUDE_DIR and MKL_LIBRARY_DIR below, 
    # e.g. on Macos uncomment
    # set(MKL_INCLUDE_DIR /Library/Frameworks/Python.framework/Versions/3.10/include)
    # set(MKL_LIBRARY_DIR /Library/Frameworks/Python.framework/Versions/3.10/lib)
    # e.g. on Windows uncomment
    # set(MKL_INCLUDE_DIR C:/Program\ Files/python310/Library/include)  
    # set(MKL_LIBRARY_DIR C:/Program\ Files/python310/Library/lib)  

    if(MKL_INCLUDE_DIR)
        message("using MKL, MKL_INCLUDE_DIR " ${MKL_INCLUDE_DIR})
    else(MKL_INCLUDE_DIR)
        message("specify MKL include directory, edit .vscode/settings.json or provide in command line mode")
        message("e.g. on Macos: -DMKL_INCLUDE_DIR=/Library/Frameworks/Python.framework/Versions/3.10/include")
        message("e.g. on Windows: -DMKL_INCLUDE_DIR=C:/Program\\ Files/python310/Library/include")
        message(FATAL_ERROR "USE_MKL=on but no MKL_INCLUDE_DIR provided")
    endif(MKL_INCLUDE_DIR)

    if(MKL_LIBRARY_DIR)
        message("using MKL, MKL_LIBRARY_DIR " ${MKL_LIBRARY_DIR})
    else(MKL_LIBRARY_DIR)
        message("specify MKL library directory, edit .vscode/settings.json or provide in command line mode")
        message("e.g. on Macos: -DMKL_LIBRARY_DIR=/Library/Frameworks/Python.framework/Versions/3.10/lib")
        message("e.g. on Windows: -DMKL_LIBRARY_DIR=C:/Program\\ Files/python310/Library/lib")
        message(FATAL_ERROR "USE_MKL=on but no MKL_LIBRARY_DIR provided")
        endif(MKL_LIBRARY_DIR)
endif (USE_MKL)


if (USE_LAPACK)
    enable_language(CXX)
    find_package(LAPACK REQUIRED)
    # manually setting alternative lapack libraries, e.g. on MACOS 
    # set(LAPACK_LIBRARIES /usr/local/opt/lapack/lib/liblapack.dylib)
    # automatically found Accelerate LAPACK is recommended!
    if(LAPACK_LIBRARIES)
    message("found LAPACK_LIBRARIES ${LAPACK_LIBRARIES}")
    else(LAPACK_LIBRARIES)
    message("did not find LAPACK, specify LAPACK_LIBRARIES in CMakeLists.txt or via -DLAPACK_LIBRARIES=<path>")
    endif(LAPACK_LIBRARIES)
endif (USE_LAPACK)


################################
# set the project name
project(GaussianElimination)

# add the executable
add_executable(GaussianElimination gauss.cpp)

target_include_directories(GaussianElimination PUBLIC ${SCLIB})
################################

################################
# set the project name
project(GaussianEff)

# add the executable
add_executable(GaussianEff gauss_efficiency.cpp)

target_include_directories(GaussianEff PUBLIC ${SCLIB})
################################
################################
# set the project name
project(GaussianStability)

# add the executable
add_executable(GaussianStability gauss_stability.cpp)

target_include_directories(GaussianStability PUBLIC ${SCLIB})
################################

################################
# set the project name
project(LUandQRFactorization)

# add the executable
add_executable(LUandQRFactorization LU_QR.cpp)

target_include_directories(LUandQRFactorization PUBLIC ${SCLIB})
################################

################################
# set the project name
project(LUandQRFactorizationC)

# add the executable
add_executable(LUandQRFactorizationC LU_QR_complex.cpp)

target_include_directories(LUandQRFactorizationC PUBLIC ${SCLIB})
################################


################################
# set the project name
project(Overdetermined)

# add the executable
add_executable(Overdetermined overdetermined.cpp)

target_include_directories(Overdetermined PUBLIC ${SCLIB})
################################

################################
# set the project name
project(SVD)

# add the executable
add_executable(SVD SVD.cpp)

target_include_directories(SVD PUBLIC ${SCLIB} ${EIGENLIB})

#if using intel mkl, link libraries (static linking, sequential, in order to generate no further dependencies)
if (USE_MKL)
    target_compile_definitions(SVD PUBLIC USE_MKL)
    target_link_directories(SVD PUBLIC ${MKL_LIBRARY_DIR})
    target_include_directories(SVD PUBLIC ${MKL_INCLUDE_DIR})
    if (UNIX)
    if (APPLE)
    target_link_libraries(SVD libmkl_intel_lp64.a libmkl_sequential.a libmkl_core.a)
    else (APPLE) ## UNIX and not APPLE is meant to run on LINUX, line below is ess. copied from MKL link instructor
    target_link_libraries(SVD  "-Wl,--start-group ${MKL_LIBRARY_DIR}/libmkl_intel_lp64.a ${MKL_LIBRARY_DIR}/libmkl_sequential.a ${MKL_LIBRARY_DIR}/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
    endif (APPLE)
    endif (UNIX)
    if (WIN32)
    target_link_libraries(SVD mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib)
    endif(WIN32)
endif (USE_MKL)
if (USE_LAPACK)
    target_compile_definitions(SVD PUBLIC USE_LAPACK)
    target_link_libraries(SVD ${LAPACK_LIBRARIES})
endif (USE_LAPACK)
################################
