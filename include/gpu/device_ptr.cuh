
#include <cstdint>
#include <cuda_runtime.h>
#include <sys/types.h>

#include "../dtypes.h"

#ifndef DEVICE_PTR__H
#define DEVICE_PTR__H

namespace tomocam::gpu {
    template <class T>
    class DevicePtr {

      private:
        T *dev_ptr_;
        dims_t dims_;

      public:
        explicit DevicePtr(dims_t dims, T *ptr) : dims_(dims), dev_ptr_(ptr) {}

        __host__ __device__ [[nodiscard]] auto dims() const { return dims_; }
        __host__ __device__ [[nodiscard]] uint_fast64_t size() const {
            return dims_.size();
        }

        // obj indexing
        __host__ __device__ T &operator[](dims_t idx3) {
            auto idx = dims_.flat_idx(idx3.x(), idx3.y(), idx3.z());
            return dev_ptr_[idx];
        }

        // const obj indexing
        __host__ __device__ const T &operator[](dims_t idx3) const {
            auto idx = dims_.flat_idx(idx3.x(), idx3.y(), idx3.z());
            return dev_ptr_[idx];
        }

        // linear indexing
        __host__ __device__ T &operator[](uint64_t idx) {
            return dev_ptr_[idx];
        }

        // const linear indexing
        __host__ __device__ const T &operator[](uint64_t idx) const {
            return dev_ptr_[idx];
        }

        // three-dim indexing
        __host__ __device__ T &operator()(int i, int j, int k) {
            auto idx = dims_.flat_idx(i, j, k);
            return dev_ptr_[idx];
        }

        // const three-dim indexing
        __host__ __device__ const T &operator()(int i, int j, int k) const {
            auto idx = dims_.flat_idx(i, j, k);
            return dev_ptr_[idx];
        }
    };

} // namespace tomocam::gpu

#endif // DEVICE_PTR__H
