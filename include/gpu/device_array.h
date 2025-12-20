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

#ifndef TOMOCAM_DEV_ARRAY__H
#define TOMOCAM_DEV_ARRAY__H

#include <cuda_runtime.h>
#include <stdexcept>

#include "array.h"
#include "gpu/cuda_check.h"
#include "gpu/gpu_memory.h"

namespace tomocam::gpu {

    template <typename T>
    class DeviceArray {
      protected:
        dims_t dims_;
        size_t size_;
        memory::cuniquePtr<T> dev_ptr_;

      public:
        DeviceArray() : dims_({0, 0, 0}), size_(0), dev_ptr_(nullptr) {}

        // Allocate space
        DeviceArray(dims_t d) : dims_(d) {
            size_ = d.n1 * d.n2 * d.n3;
            dev_ptr_ = memory::make_cuniquePtr<T>(size_);
        }

        /* create device array from Slice */
        DeviceArray(const Slice<T> &rhs)
            : dims_(rhs.dims()), size_(rhs.size()),
              dev_ptr_(memory::make_cuniquePtr<T>(rhs.size())) {
            SAFE_CALL(cudaMemcpy(dev_ptr_.get(), rhs.begin(), rhs.bytes(),
                                 cudaMemcpyHostToDevice));
        }

        //  copy constructor
        DeviceArray(const DeviceArray<T> &rhs)
            : dims_(rhs.dims_), size_(rhs.size_),
              dev_ptr_(memory::make_cuniquePtr<T>(size_)) {
            SAFE_CALL(cudaMemcpy(dev_ptr_.get(), rhs.dev_ptr_.get(), rhs.bytes(),
                                 cudaMemcpyDeviceToDevice));
        }

        // assignment operator
        DeviceArray<T> &operator=(const DeviceArray &rhs) {
            if (this != &rhs) {
                auto new_ptr = memory::make_cuniquePtr<T>(rhs.size_);
                SAFE_CALL(cudaMemcpy(new_ptr.get(), rhs.dev_ptr_.get(), rhs.bytes(),
                                     cudaMemcpyDeviceToDevice));
                dims_ = rhs.dims_;
                size_ = rhs.size_;
                dev_ptr_ = std::move(new_ptr);
            }
            return *this;
        }

        // cuniquePtr makes destructor redundant
        ~DeviceArray() = default;

        // move constructor
        DeviceArray(DeviceArray<T> &&rhs) = default;

        // move assignment operator
        DeviceArray<T> &operator=(DeviceArray<T> &&rhs) = default;

        // access to the beginning of the array
        T *begin() { return dev_ptr_.get(); };
        const T *begin() const { return dev_ptr_.get(); };

        // access to the end of the array
        T *end() { return dev_ptr_.get() + size_; };
        const T *end() const { return dev_ptr_.get() + size_; };

        // size of the array
        [[nodiscard]] size_t size() const { return size_; }

        // bytes of the array
        [[nodiscard]] size_t bytes() const { return sizeof(T) * size_; }

        // get array dims
        [[nodiscard]] dims_t dims() const { return dims_; }

        // get number of slices
        [[nodiscard]] size_t nslices() const { return dims_.n1; }

        // get number of rows
        [[nodiscard]] size_t nrows() const { return dims_.n2; }

        // get number of columns
        [[nodiscard]] size_t ncols() const { return dims_.n3; }
    };

} // namespace tomocam::gpu

#endif // TOMOCAM_DEV_ARRAY__H
