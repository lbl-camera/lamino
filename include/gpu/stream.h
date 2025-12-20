
#ifndef TOMOCAM_RAII_UTILS__H
#define TOMOCAM_RAII_UTILS__H

#include <cuda_runtime.h>

#include "gpu_memory.h"

namespace tomocam::gpu {
    // RAII wrapper for CUDA stream
    class Stream {
      private:
        cudaStream_t stream_;

      public:
        Stream() { cudaStreamCreate(&stream_); }
        ~Stream() { cudaStreamDestroy(stream_); }
        cudaStream_t get() const { return stream_; }
    };

    // RAII wrapper for pinned host memory
    class PinnedBuffer {
      private:
        memory::pinned_ptr<void> ptr_;
        size_t size_;

      public:
        PinnedBuffer(size_t size) : size_(size) {
            ptr_ = memory::make_pinned_ptr<void>(size);
        }
        void *get() const { return ptr_.get(); }
        size_t size() const { return size_; }
    };
} // namespace tomocam::gpu
#endif // TOMOCAM_RAII_UTILS__H
