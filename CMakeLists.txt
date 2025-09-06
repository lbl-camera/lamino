
cmake_minimum_required(VERSION 3.12)

project(Testbed CXX)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(SRC
    src/projection.cpp
    src/results.cpp
    src/fftconv.cpp
    src/mbir1.cpp
   )
include_directories(${CMAKE_SOURCE_DIR}/src)

camke_export_compile_commands(ON)

# C++ std
set(CMAKE_CXX_STANDARD 23)

find_package(OpenMP)
if(OpenMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OpenMP_FOUND)

# tbb
find_package(TBB REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)

# finufft
find_package(finufft REQUIRED)

add_library(tomocam SHARED ${SRC})
target_link_libraries(tomocam PUBLIC fftw3_omp)
target_link_libraries(tomocam PUBLIC finufft::finufft)
target_link_libraries(tomocam PUBLIC TBB::tbb)


option(ENABLE_TESTS "Enable tests" ON)
if (${ENABLE_TESTS})
    add_subdirectory(test)
endif()
