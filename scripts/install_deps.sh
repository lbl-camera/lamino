#!/bin/bash
# System package dependencies for building XMCD/lamino (tomocam).
#
# Target OS: Ubuntu 24.04 (Noble)
# Assumes CUDA Toolkit is already installed separately (tested with 13.1,
# GPU compute capability 8.9 / Ada Lovelace).
#
# NOT covered here (built from source, not apt packages):
#   - FINUFFT / cuFINUFFT — see https://github.com/flatironinstitute/finufft
#     Build into ~/finufft (see cmake/Findfinufft.cmake for the search path),
#     e.g.:
#       git clone https://github.com/flatironinstitute/finufft.git ~/finufft
#       cd ~/finufft
#       cmake -S . -B build -GNinja \
#           -DCMAKE_INSTALL_PREFIX=$HOME/finufft \
#           -DCMAKE_INSTALL_LIBDIR=lib64 \
#           -DCMAKE_CXX_COMPILER=clang++ \
#           -DCMAKE_C_COMPILER=clang \
#           -DCMAKE_CUDA_COMPILER=nvcc \
#           -DCMAKE_CUDA_HOST_COMPILER=clang \
#           -DCMAKE_BUILD_TYPE=Release \
#           -DBUILD_TESTS:BOOL=OFF -DBUILD_EXAMPLES:BOOL=OFF \
#           -DFINUFFT_USE_CUDA:BOOL=ON \
#           -DCMAKE_CUDA_ARCHITECTURES=89 \
#           -DFINUFFT_STATIC_LINKING:BOOL=OFF \
#           -DFINUFFT_BUILD_PYTHON:BOOL=OFF
#       cmake --build build && cmake --install build
#
# NOTE: CMAKE_CUDA_ARCHITECTURES in the top-level CMakeLists.txt is hardcoded
# to "60;61;70;75;80;86", which includes Pascal/Volta archs (60/61/70) that
# CUDA 13.x can no longer compile for. Update it to match the target GPU
# (e.g. "89" for compute capability 8.9) before building.

set -euo pipefail

sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    ninja-build \
    git \
    clang-18 \
    libomp-18-dev \
    libtbb-dev \
    libfftw3-dev \
    libtiff-dev
