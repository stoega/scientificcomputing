cmake_minimum_required(VERSION 3.13)
set (CMAKE_CXX_STANDARD 14)

set(SCLIB ../../SClinalg/)
set(EIGENLIB ../../eigen/)

option( USE_MKL "enable MKL" off)


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
else (USE_MKL)
    enable_language(CXX)
    find_package(LAPACK QUIET)
    # manually setting alternative lapack libraries, e.g. on MACOS 
    # set(LAPACK_LIBRARIES /usr/local/opt/lapack/lib/liblapack.dylib)
    # automatically found Accelerate LAPACK is recommended!
    if(LAPACK_LIBRARIES)
    message("found LAPACK_LIBRARIES ${LAPACK_LIBRARIES}")
    else(LAPACK_LIBRARIES)
    message("did not find LAPACK, specify LAPACK_LIBRARIES in .vscode/settings.json or CMakeLists.txt or via -DLAPACK_LIBRARIES=<path>")
    endif(LAPACK_LIBRARIES)
endif (USE_MKL)

################################
# set the project name
project(LAPACKTEST)

# add the executable
add_executable(LAPACKTEST LAPACK.cpp)

target_include_directories(LAPACKTEST PUBLIC ${SCLIB} ${EIGENLIB})

if (USE_MKL)
    target_compile_definitions(LAPACKTEST PUBLIC USE_MKL)
    target_link_directories(LAPACKTEST PUBLIC ${MKL_LIBRARY_DIR})
    target_include_directories(LAPACKTEST PUBLIC ${MKL_INCLUDE_DIR})
    if (UNIX)
    if (APPLE)
    # static linking, sequential
    target_link_libraries(LAPACKTEST libmkl_intel_lp64.a libmkl_sequential.a libmkl_core.a)
    # dynamic linking single lib
    #target_link_libraries(LAPACKTEST libmkl_rt.2.dylib)
    else (APPLE) ## UNIX and not APPLE means LINUX here
    # static linking, sequential
    target_link_libraries(LAPACKTEST  "-Wl,--start-group ${MKL_LIBRARY_DIR}/libmkl_intel_lp64.a ${MKL_LIBRARY_DIR}/libmkl_sequential.a ${MKL_LIBRARY_DIR}/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")
    endif(APPLE)
    endif (UNIX)
    if (WIN32)
    # static linking, sequential
    target_link_libraries(LAPACKTEST mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib)
    # dynamic linking single library
    # target_link_libraries(LAPACKTEST mkl_rt.lib)
    endif(WIN32)
else (USE_MKL)
    target_link_libraries(LAPACKTEST "${LAPACK_LIBRARIES}")
endif (USE_MKL)
################################
