
cmake_minimum_required(VERSION 3.8)

project(Testbed CXX)

set(SRC
    fft.cpp
    deapod.cpp
    convolve.cpp
    projection.cpp
    results.cpp
    main.cpp
   )

# C++11 Stuff
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP)
if(OpenMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OpenMP_FOUND)

add_executable(nufft ${SRC})
target_link_libraries(nufft fftw3f)
