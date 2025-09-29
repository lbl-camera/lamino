
cmake_minimum_required(VERSION 3.20)

project(Testbed CXX)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Add cuda
option(USE_GPU "Use GPU Acceleration" TRUE)
if (${USE_GPU})
  enable_language(CUDA)
endif()

# C++ std
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC
    src/projection.cpp
    src/results.cpp
   )
if (${USE_GPU})
  list(APPEND SRC src/gpu/polar_grid.cu)
endif()
include_directories(${CMAKE_SOURCE_DIR}/include)

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
    enable_testing()
    add_subdirectory(tests)
endif()
