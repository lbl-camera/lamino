
cmake_minimum_required(VERSION 3.20)

project(Testbed CXX)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# C++ std
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC
    src/projection.cpp
    src/results.cpp
   )
include_directories(${CMAKE_SOURCE_DIR}/src)

find_package(OpenMP)
if(OpenMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OpenMP_FOUND)

# tbb
find_package(TBB REQUIRED)
find_package(FFTW REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)

# finufft
find_package(finufft REQUIRED)
add_library(tomocam SHARED ${SRC})
target_link_libraries(tomocam PUBLIC FFTW::double)
target_link_libraries(tomocam PUBLIC FFTW::float)
target_link_libraries(tomocam PUBLIC finufft::finufft)
target_link_libraries(tomocam PUBLIC TBB::tbb)
set_property(TARGET tomocam PROPERTY CXX_STANDARD 23)


option(ENABLE_TESTS "Enable tests" OFF)
if (${ENABLE_TESTS})
    add_subdirectory(test)
endif()
