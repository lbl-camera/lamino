#! /bin/bash

cmake -S . -B gpu -G Ninja  \
    -DCMAKE_BUILD_TYPE=Release  \
    -DCMAKE_CXX_COMPILER=clang++ \
    -DCMAKE_CUDA_HOST_COMPILER=clang++ \
    -DENABLE_TESTS:BOOL=ON \
    -DUSE_CUDA:BOOL=ON \
    -DCMAKE_CXX_FLAGS="-fno-omit-frame-pointer"

cmake --build gpu --target all -j 16
