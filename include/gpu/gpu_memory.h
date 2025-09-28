#include <iostream>
#include <memory>
#include <stdexcept>

#include <cuda_runtime.h>
#include <string>

#ifndef GPU_MEMORY__H
#define GPU_MEMORY__H

namespace tomocam::gpu {
    struct cudaDelete {
        void operator()(void *ptr) const noexcept {
            if (ptr != nullptr) {
                auto err = cudaFree(ptr);
                if (err != cudaSuccess) {
                    std::cerr << "cudaFree failed: " << cudaGetErrorString(err)
                              << '\n';
                }
            }
        }
    };

    template <typename T>
    using cudaPtr = std::unique_ptr<T, cudaDelete>;

    template <typename T>
    cudaPtr<T> make_cudaPtr(size_t count) {
        T *raw = nullptr;
        auto err = cudaMallocManaged(&raw, sizeof(T) * count);
        if (err != cudaSuccess) {
            throw std::runtime_error(
                std::string("failed to allocated unified memory"));
        }
        return cudaPtr<T>(raw);
    }
} // namespace tomocam::gpu

#endif // GPU_MEMORY__H
