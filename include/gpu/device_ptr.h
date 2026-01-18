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

#ifndef DEVICE_PTR_H
#define DEVICE_PTR_H

#include <cstdint>
#include <cuda_runtime.h>

#include "dtypes.h"

namespace tomocam::gpu {

    /// Concept for containers compatible with DevicePtr wrapping.
    /// Requires containers to provide:
    /// - begin() method returning a pointer to value_type
    /// - dims() method returning dims_t with dimensions information
    template <class C>
    concept WappableToDevicePtr = requires(C c) {
        { c.begin() } -> std::convertible_to<typename C::value_type *>;
        { c.dims() } -> std::convertible_to<dims_t>;
    };

    /// Non-owning view wrapper for GPU device pointers with 3D indexing support.
    /// Designed for implicit construction in CUDA kernel invocations. Provides
    /// convenient multi-dimensional indexing for GPU data while maintaining a
    /// lightweight wrapper around raw device pointers.
    ///
    /// @tparam T Element type of the device array
    ///
    /// Usage:
    /// @code
    /// __global__ void my_kernel(DevicePtr<float> data) {
    ///     float val = data(i, j, k);  // 3D indexing
    ///     float val2 = data[i*n2*n3 + j*n3 + k];  // Linear indexing
    /// }
    ///
    /// DeviceArray<float> dev_array(...);
    /// my_kernel<<<blocks, threads>>>(dev_array);  // Implicit conversion
    /// @endcode
    template <typename T>
    class DevicePtr {
      private:
        T *dev_ptr_;
        dims_t dims_;

        __host__ __device__ size_t flat_idx(int i, int j, int k) const {
            return (static_cast<size_t>(i) * dims_.n2 * dims_.n3) +
                   (static_cast<size_t>(j) * dims_.n3) + static_cast<size_t>(k);
        }

      public:
        /// Implicitly constructs DevicePtr from compatible containers in kernel
        /// calls.
        template <WappableToDevicePtr C>
        __host__ DevicePtr(C &container)
            : dev_ptr_(container.begin()), dims_(container.dims()) {}

        /// Returns dimensions of the underlying 3D array.
        __host__ __device__ [[nodiscard]] auto dims() const { return dims_; }

        /// Returns total number of elements in the array.
        __host__ __device__ [[nodiscard]] size_t size() const {
            return (dims_.n1 * dims_.n2 * dims_.n3);
        }

        /// Accesses element via 3D index struct.
        __host__ __device__ T &operator[](uint3 idx3) {
            auto idx = flat_idx(idx3.x, idx3.y, idx3.z);
            return dev_ptr_[idx];
        }

        /// Accesses element via 3D index struct (const version).
        __host__ __device__ const T &operator[](uint3 idx3) const {
            auto idx = flat_idx(idx3.x, idx3.y, idx3.z);
            return dev_ptr_[idx];
        }

        /// Accesses element via linear index.
        __host__ __device__ T &operator[](size_t idx) { return dev_ptr_[idx]; }

        /// Accesses element via linear index (const version).
        __host__ __device__ const T &operator[](size_t idx) const {
            return dev_ptr_[idx];
        }

        /// Accesses element via 3D indices (i, j, k).
        __host__ __device__ T &operator()(int i, int j, int k) {
            auto idx = flat_idx(i, j, k);
            return dev_ptr_[idx];
        }

        /// Accesses element via 3D indices (i, j, k) (const version).
        __host__ __device__ const T &operator()(int i, int j, int k) const {
            auto idx = flat_idx(i, j, k);
            return dev_ptr_[idx];
        }
    };

} // namespace tomocam::gpu

#endif // DEVICE_PTR_H
