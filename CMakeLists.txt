
cmake_minimum_required(VERSION 3.8)

project(Testbed CXX)

set(SRC
    src/fft.cpp
    src/deapod.cpp
    src/fftconv.cpp
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

# flat iron nufft
set(FINUFFT_DIR "/home/dkumar/finufft")
include_directories(${FINUFFT_DIR}/include)
link_directories(${FINUFFT_DIR}/lib)

add_library(nufft STATIC ${SRC})
target_link_libraries(nufft PRIVATE fftw3f)

add_executable(backward test/backward.cpp)
target_link_libraries(backward PRIVATE nufft)
target_link_libraries(backward PRIVATE finufft)

add_executable(forward test/forward.cpp)
target_link_libraries(forward PRIVATE nufft)
target_link_libraries(forward PRIVATE finufft)

add_executable(sirt test/main.cpp)
target_link_libraries(sirt PRIVATE nufft)
target_link_libraries(sirt PRIVATE finufft)

add_executable(ffterr test/transform_error.cpp)
target_link_libraries(ffterr PRIVATE nufft)
target_link_libraries(ffterr PRIVATE finufft)

add_executable(test_array test/array.cpp)
