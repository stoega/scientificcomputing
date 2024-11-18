cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 11)

## avoid linker error when using XCode15, see https://developer.apple.com/forums/thread/737707
# if (APPLE)
# set(CMAKE_CXX_FLAGS "-no-pie")
# endif (APPLE)

################################
# set the project name
project(BasicTypes)

# add the executable
add_executable(BasicTypes basictypes.cpp)

################################


################################
# set the project name
project(Functions)

# add the executable
add_executable(Functions functions.cpp)

################################


################################
# set the project name
project(RefsPtrs)

# add the executable
add_executable(RefsPtrs references_pointers.cpp)

################################

################################
# set the project name
project(CallByRef)

# add the executable
add_executable(CallByRef cbr.cpp)

################################

################################
# set the project name
project(Memory)

# add the executable
add_executable(Memory memory.cpp)

################################

################################
# set the project name
project(Namespaces)

# add the executable
add_executable(Namespaces namespaces.cpp)

################################


################################
# set the project name
project(IOStream)

# add the executable
add_executable(IOStream iostream.cpp)

################################

################################
# set the project name
project(Complex)

# add the executable
add_executable(Complex complex.cpp)

################################
