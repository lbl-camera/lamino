
cmake_minimum_required(VERSION 3.8)

project(Testbed CXX)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(SRC
    src/fft.cpp
    src/projection.cpp
    src/proj2.cpp
    src/results.cpp
    src/fftconv.cpp
   )
include_directories(${CMAKE_SOURCE_DIR}/src)

# C++17 Stuff
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

find_package(OpenMP)
if(OpenMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OpenMP_FOUND)

find_package(HDF5 REQUIRED COMPONENTS CXX)

# finufft
find_package(finufft REQUIRED)

add_library(tomocam SHARED ${SRC})
target_link_libraries(tomocam PRIVATE fftw3)
target_link_libraries(tomocam PRIVATE finufft::finufft)

add_executable(sirt test/main.cpp)
target_link_libraries(sirt PRIVATE ${HDF5_LIBRARIES})
target_link_libraries(sirt PRIVATE tomocam)
target_link_libraries(sirt PRIVATE finufft::finufft)

#add_executable(toeplitz test/toeplitz.cpp)
#target_link_libraries(toeplitz PRIVATE ${HDF5_LIBRARIES})
#target_link_libraries(toeplitz PRIVATE tomocam)
#target_link_libraries(toeplitz PRIVATE finufft::finufft)
