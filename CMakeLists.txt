
cmake_minimum_required(VERSION 3.20)

project(Testbed CXX)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# C++ std
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DDEBUG)
endif()
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -march=native")
endif()

set(SRC
    src/bregman.cpp
    src/finitediff.cpp
    src/conjgrad.cpp
    src/projection.cpp
    src/gradient.cpp
    src/mbir_sb.cpp
    src/mbir_cg.cpp
    src/nagopt.cpp
    src/qggmrf.cpp
    src/recon.cpp
   )
include_directories(${CMAKE_SOURCE_DIR}/include)

find_package(OpenMP)
if(OpenMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OpenMP_FOUND)

# toml++
include(FetchContent)
FetchContent_Declare(
  tomlplusplus
  GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(tomlplusplus)
include_directories(${tomlplusplus_SOURCE_DIR}/include)


# tbb
find_package(TBB REQUIRED)
find_package(FFTW REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(TIFF REQUIRED)

# finufft
find_package(finufft REQUIRED)
add_library(tomocam SHARED ${SRC})
target_link_libraries(tomocam PUBLIC FFTW::double)
target_link_libraries(tomocam PUBLIC FFTW::float)
target_link_libraries(tomocam PUBLIC finufft::finufft)
target_link_libraries(tomocam PUBLIC TBB::tbb)
target_link_libraries(tomocam PUBLIC TIFF::TIFF)
set_property(TARGET tomocam PROPERTY CXX_STANDARD 23)

# recon executable
add_executable(recon src/recon.cpp)
target_link_libraries(recon PRIVATE tomocam)

# install executables and library
install(TARGETS recon DESTINATION bin)
install(TARGETS tomocam DESTINATION lib)

option(ENABLE_TESTS "Enable tests" OFF)
if (${ENABLE_TESTS})
    enable_testing()
    add_subdirectory(tests)
endif()
