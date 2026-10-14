#! /bin/bash

cmake -S . -B cpu_debug -G Ninja  \
    -DCMAKE_BUILD_TYPE=Debug  \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DENABLE_TESTS:BOOL=ON \
    -DUSE_CUDA:BOOL=OFF \
    -DCMAKE_CXX_FLAGS="-fno-omit-frame-pointer" \
    -DOpenMP_C_FLAGS="-fopenmp" \
    -DOpenMP_CXX_FLAGS="-fopenmp" \
    -DOpenMP_C_LIB_NAMES="omp" \
    -DOpenMP_CXX_LIB_NAMES="omp" \
    -DOpenMP_omp_LIBRARY=/usr/lib/llvm-18/lib/libomp.so

cmake --build cpu_debug --target all -j 16
