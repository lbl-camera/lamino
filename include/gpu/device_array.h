/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 *National Laboratory (subject to receipt of any required approvals from the
 *U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 *IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 *the U.S. Government has been granted for itself and others acting on its
 *behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 *to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#ifndef TOMOCAM_DEV_ARRAY_H
#define TOMOCAM_DEV_ARRAY_H

#include <cuda_runtime.h>
#include <stdexcept>

#include "array.h"
#include "gpu/gpu_memory.h"
#include "gpu/utils.h"

namespace tomocam::gpu {

    /// Container wrapping a CUDA device pointer as a 3D array.
    /// Manages device memory with RAII semantics using cunique_ptr. Provides
    /// 3D indexing through implicit conversion to DevicePtr for use in CUDA kernels.
    ///
    /// @tparam T Element type of the device array
    ///
    /// Usage:
    /// @code
    /// DeviceArray<float> dev_array({10, 20, 30});  // Allocate 10x20x30
    ///
    /// Slice<float> host_data(...);
    /// DeviceArray<float> dev_from_host(host_data);  // Copy from host
    ///
    /// __global__ void my_kernel(DevicePtr<float> data) {
    ///     float val = data(i, j, k);
    /// }
    /// my_kernel<<<blocks, threads>>>(dev_array);  // Implicit conversion to
    /// DevicePtr
    /// @endcode
    template <typename T>
    class DeviceArray {
      public:
        using value_type = T;

      protected:
        dims_t dims_;
        size_t size_;
        memory::cunique_ptr<T> dev_ptr_;

      public:
        /// Default constructor creates empty array.
        DeviceArray() : dims_({0, 0, 0}), size_(0), dev_ptr_(nullptr) {}

        /// Allocates device memory for a 3D array of given dimensions.
        DeviceArray(dims_t d) : dims_(d) {
            size_ = d.n1 * d.n2 * d.n3;
            dev_ptr_ = memory::make_cunique_ptr<T>(size_);
        }

        /// Constructs device array from host Slice, copying data to device.
        DeviceArray(const Slice<T> &rhs)
            : dims_(rhs.dims()), size_(rhs.size()),
              dev_ptr_(memory::make_cunique_ptr<T>(rhs.size())) {
            SAFE_CALL(cudaMemcpy(dev_ptr_.get(), rhs.begin(), rhs.bytes(),
                                 cudaMemcpyHostToDevice));
        }

        /// Copy constructor is deleted to prevent unintended device memory
        /// duplication.
        DeviceArray(const DeviceArray<T> &rhs) = delete;
        /// Copy assignment is deleted to prevent unintended device memory
        /// duplication.
        DeviceArray<T> &operator=(const Slice<T> &rhs) = delete;

        /// Destructor automatically managed by cunique_ptr.
        ~DeviceArray() = default;

        /// Creates explicit clone of the device array on GPU.
        DeviceArray<T> clone() const {
            DeviceArray<T> copy(dims_);
            SAFE_CALL(cudaMemcpy(copy.dev_ptr_.get(), dev_ptr_.get(), bytes(),
                                 cudaMemcpyDeviceToDevice));
            return copy;
        }

        /// Move constructor for efficient transfer of ownership.
        DeviceArray(DeviceArray<T> &&rhs) = default;

        /// Move assignment operator for efficient transfer of ownership.
        DeviceArray<T> &operator=(DeviceArray<T> &&rhs) = default;

        /// Returns pointer to beginning of device memory.
        T *begin() { return dev_ptr_.get(); };
        /// Returns const pointer to beginning of device memory.
        const T *begin() const { return dev_ptr_.get(); };

        /// Returns pointer to one past the end of device memory.
        T *end() { return dev_ptr_.get() + size_; };
        /// Returns const pointer to one past the end of device memory.
        const T *end() const { return dev_ptr_.get() + size_; };

        /// Returns total number of elements in the array.
        [[nodiscard]] size_t size() const { return size_; }

        /// Returns total size in bytes of the array.
        [[nodiscard]] size_t bytes() const { return sizeof(T) * size_; }

        /// Returns dimensions (n1, n2, n3) of the 3D array.
        [[nodiscard]] dims_t dims() const { return dims_; }

        /// Returns number of slices (first dimension).
        [[nodiscard]] size_t nslices() const { return dims_.n1; }

        /// Returns number of rows (second dimension).
        [[nodiscard]] size_t nrows() const { return dims_.n2; }

        /// Returns number of columns (third dimension).
        [[nodiscard]] size_t ncols() const { return dims_.n3; }
    };

} // namespace tomocam::gpu

#endif // TOMOCAM_DEV_ARRAY__H
