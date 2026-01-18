
#ifndef TOMOCAM_GPU_UTILS_H
#define TOMOCAM_GPU_UTILS_H

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#include "dtypes.h"

#define SAFE_CALL(call)                                                             \
    {                                                                               \
        cudaError_t err = call;                                                     \
        if (err != cudaSuccess) {                                                   \
            fprintf(stderr, "CUDA error in file '%s' in line %i : %s.\n", __FILE__, \
                    __LINE__, cudaGetErrorString(err));                             \
            exit(EXIT_FAILURE);                                                     \
        }                                                                           \
    }

#ifdef __CUDACC__
__device__ inline bool operator<(const uint3 &a, const tomocam::dims_t &b) {
    if (b > a) {
        return true;
    } else {
        return false;
    }
}

inline dim3 make_grid(const tomocam::dims_t &dims, const dim3 &threads) {
    return dim3((dims.n1 + threads.x - 1) / threads.x,
                (dims.n2 + threads.y - 1) / threads.y,
                (dims.n3 + threads.z - 1) / threads.z);
}

__device__ inline uint3 Index3D() {
    uint3 idx;
    idx.x = blockIdx.x * blockDim.x + threadIdx.x;
    idx.y = blockIdx.y * blockDim.y + threadIdx.y;
    idx.z = blockIdx.z * blockDim.z + threadIdx.z;
    return idx;
}
#endif

#endif // TOMOCAM_GPU_UTILS_H
