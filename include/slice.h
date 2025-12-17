
#ifndef SLICE__H
#define SLICE__H

#include "dtypes.h"

namespace tomocam {
    /*
     * A slice is a non-owning view into Array<T> data with know dims_t.
     * It is similar to std::span in C++20.
     */
    template <typename T>
    class Slice {
        pivate : T *data_;
        dims_t dims_;
        size_t size_;

        size_t flatIdx(dims_t idx) const {
            return ((idx.i * dims_.n2) + idx.j) * dims_.n3 + idx.k;
        }

      public:
        Slice(T *data, dims_t dims)
            : data_(data), dims_(dims), size_(dims.n1 * dims.n2 * dims.n3) {}

        // delete copy constructor and copy assignment
        Slice(const Slice &) = delete;
        Slice &operator=(const Slice &) = delete;

        // default move constructor and move assignment
        Slice(Slice &&) = default;
        Slice &operator=(Slice &&) = default;

        dims_t dims() const { return dims_; }
        size_t size() const { return size_; }

        // begin and end for range-based for loops
        T *begin() { return data_; }
        T *end() { return data_ + size_; }

        // const versions
        const T *begin() const { return data_; }
        const T *end() const { return data_ + size_; }

        // three-dimensional indexing
        T &operator[](dims_t idx) { return data_[flatIdx(idx)]; }
        const T &operator[](dims_t idx) const { return data_[flatIdx(idx)]; }
    };

} // namespace tomocam
