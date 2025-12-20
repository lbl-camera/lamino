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
#ifndef GPUMEMORY__H
#define GPUMEMORY__H

#include <iostream>
#include <memory>
#include <stdexcept>

#include <cuda_runtime.h>
#include <string>

namespace tomocam::gpu::memory {
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
    using cunique_ptr = std::unique_ptr<T, cudaDelete>;

    template <typename T>
    cunique_ptr make_cunique_ptr(std::size_t count) {
        T *raw = nullptr;
        auto err = cudaMalloc(&raw, sizeof(T) * count);
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("failed to allocated gpu memory"));
        }
        return cunique_ptr<T>(raw);
    }

    struct pinnedDelete {
        void operator()(void *ptr) const noexcept {
            if (ptr != nullptr) {
                auto err = cudaFreeHost(ptr);
                if (err != cudaSuccess) {
                    std::cerr << "cudaFreeHost failed: " << cudaGetErrorString(err)
                              << '\n';
                }
            }
        }
    };

    template <typename T>
    using pinned_ptr = std::unique_ptr<T, pinnedDelete>;

    template <typename T>
    pinned_ptr<T> make_pinned_ptr(std::size_t count) {
        T *raw = nullptr;
        auto err = cudaMallocHost(&raw, sizeof(T) * count);
        if (err != cudaSuccess) {
            throw std::runtime_error(
                std::string("failed to allocated pinned host memory"));
        }
        return pinned_ptr<T>(raw);
    }
} // namespace tomocam::gpu::memory

#endif // GPUMEMORY__H
