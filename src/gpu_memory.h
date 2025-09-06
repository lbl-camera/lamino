#include <iostream>
#include <memory>
#include <stdexcept>

#include <cuda_runtime.h>
#include <string>

#ifndef GPU_MEMORY__H
#define GPU_MEMORY__H

namespace tomocam {
    namespace gpu {
        struct gpuDeleter {
            void operator()(void *ptr) const noexcept {
                if (ptr) {
                    auto err = cudaFree(ptr);
                    if (err != cudaSuccess) {
                        std::cerr
                            << "cudaFree failed: " << cudaGetErrorString(err)
                            << std::endl;
                    }
                }
            }
        };

        template <typename T>
        using uniquePtr = std::unique_ptr<T, gpuDeleter>;

        template <typename T>
        uniquePtr<T> make_uniquePtr(size_t count) {
            T *raw = nullptr;
            auto err = cudaMallocManaged(&raw, sizeof(T) * count);
            if (err != cudaSuccess) {
                throw std::runtime_error(
                    std::string("failed to allocated unified memory"));
            }
            return uniquePtr<T>(raw);
        }
    } // namespace gpu
} // namespace tomocam

#endif // GPU_MEMORY__H
