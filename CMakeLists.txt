
cmake_minimum_required(VERSION 3.8)

project(Testbed CXX)

set(SRC
    src/fft.cpp
    src/deapod.cpp
    src/convolve.cpp
    src/projection.cpp
    src/results.cpp
   )
include_directories(${CMAKE_SOURCE_DIR}/src)

# C++17 Stuff
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

find_package(OpenMP)
if(OpenMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OpenMP_FOUND)


add_library(nufft STATIC ${SRC})
target_link_libraries(nufft PRIVATE fftw3f)

add_executable(sirt test/main.cpp)
target_link_libraries(sirt nufft)

add_executable(gridrec test/gridrec.cpp)
target_link_libraries(gridrec nufft)
