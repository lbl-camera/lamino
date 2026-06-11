
cmake_minimum_required(VERSION 3.20)

project(XMCD CXX)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
option("USE_CUDA" "Use GPU acceleration" OFF)
message(STATUS "USE_CUDA: ${USE_CUDA}")
if (NOT USE_CUDA)
    message(STATUS "use -DUSE_CUDA:BOOL=ON to enable GPU acceleration")
endif()
if (USE_CUDA)
    enable_language(CUDA)
endif()


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
   )
include_directories(${CMAKE_SOURCE_DIR}/include)

# append cu files if CUDA is enabled
if (USE_CUDA)
    list(APPEND SRC
        src/gpu/mbir.cu
        src/gpu/device_array.cu
        src/gpu/projection.cu
        src/gpu/gradient.cu
        src/gpu/polar_grid.cu
        src/gpu/array_ops.cu
        src/gpu/bregman.cu
        src/gpu/conjgrad.cu
        src/gpu/fftshift.cu
        src/gpu/finitediff.cu
        src/gpu/padding.cu
        src/gpu/xmcdprojs.cu
    )
    include_directories(${CMAKE_SOURCE_DIR}/src/cuda)
endif()

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

if (USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    set(CMAKE_CUDA_ARCHITECTURES "60;61;70;75;80;86")
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")
    add_definitions(-DUSE_CUDA)
endif()


# tbb
find_package(TBB REQUIRED)
find_package(FFTW REQUIRED)
find_package(TIFF REQUIRED)

# finufft
find_package(finufft REQUIRED)

if (CMAKE_CUDA_COMPILER AND NOT TARGET finufft::cufinufft)
    message(FATAL_ERROR "CUDA is enabled but finufft::cufinufft target not found. Ensure cufinufft library is installed.")
endif()

add_library(tomocam SHARED ${SRC})
target_link_libraries(tomocam PUBLIC FFTW::double)
target_link_libraries(tomocam PUBLIC FFTW::float)
target_link_libraries(tomocam INTERFACE finufft::finufft)
if (TARGET finufft::cufinufft)
    target_link_libraries(tomocam INTERFACE finufft::cufinufft)
endif()
target_link_libraries(tomocam PUBLIC TBB::tbb)
target_link_libraries(tomocam PUBLIC TIFF::TIFF)
if (USE_CUDA)
    target_link_libraries(tomocam PUBLIC CUDA::cufft)
    target_link_libraries(tomocam PUBLIC CUDA::cudart)
endif()
set_property(TARGET tomocam PROPERTY CXX_STANDARD 23)

# recon executable
add_executable(recon src/recon.cpp)
target_link_libraries(recon PRIVATE tomocam)

# forward executable
add_executable(forward src/forward.cpp)
target_link_libraries(forward PUBLIC tomocam)
target_link_libraries(forward PUBLIC TIFF::TIFF)

# adjoint executable
add_executable(adjoint src/adjoint.cpp)
target_link_libraries(adjoint PUBLIC tomocam)
target_link_libraries(adjoint PUBLIC TIFF::TIFF)

# install executables and library
install(TARGETS recon forward adjoint DESTINATION bin)
install(TARGETS tomocam DESTINATION lib)

option(ENABLE_TESTS "Enable tests" OFF)
if (${ENABLE_TESTS})
    enable_testing()
    add_subdirectory(tests)
endif()
