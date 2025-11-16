// clang-format off
/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 * National Laboratory (subject to receipt of any required approvals from the
 * U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 * IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 * the U.S. Government has been granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 * to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */
 //clang-format on
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
