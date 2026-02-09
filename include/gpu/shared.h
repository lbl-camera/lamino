
#ifndef CUDA_SHARED_H
#define CUDA_SHARED_H

#include <cuda_runtime.h>

namespace tomo::gpu {
    template <typename T>
    struct SharedMemory3D {
        const uint3 dims;
        T *data;

        __device__ SharedMemory3D(uint3 dimensions = make_uint3(3, 16, 16))
            : dims(dimensions) {
            extern __shared__ __align__(sizeof(T)) unsigned char s_data[];
            data = reinterpret_cast<T *>(s_data);
        }
        __device__ T &operator[](uint3 idx) {
#ifdef DEBUG
            assert(idx.x < dims.x && idx.y < dims.y && idx.z < dims.z);
#endif
            return data[idx.x * dims.y * dims.z + idx.y * dims.z + idx.z];
        }
        __device__ const T &operator[](uint3 idx) const {
#ifdef DEBUG
            assert(idx.x < dims.x && idx.y < dims.y && idx.z < dims.z);
#endif
            return data[idx.x * dims.y * dims.z + idx.y * dims.z + idx.z];
        }
        __device__ T &operator[](size_t idx) {
#ifdef DEBUG
            assert(idx < dims.x * dims.y * dims.z);
#endif
            return data[idx];
        }
        __device__ const T &operator[](size_t idx) const {
#ifdef DEBUG
            assert(idx < dims.x * dims.y * dims.z);
#endif
            return data[idx];
        }

        __device__ size_t size() const {
            return dims.x * dims.y * dims.z * sizeof(T);
        }
    };
} // namespace tomo::gpu

#endif // CUDA_SHARED_H
