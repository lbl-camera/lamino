
#ifndef ARRAY3__H
#define ARRAY3__H

#include <string.h>

#include "dtypes.h"

template <typename T>
class Array {
    private:
        dim3_t dims_;
        size_t size_;
        T * ptr_;

    public:
        Array(int x, int y, int z) {
            dims_ = {x, y, z};
            size_ = x * y * z;
            ptr_ = new T[size_];
        }
        ~Array() {
            if (ptr_) delete [] ptr_;
        }

        Array(const Array &rhs) = delete;
        Array& operator=(const Array &rhs) = delete;

        void clear() { memset(ptr_, 0, size_ * sizeof(T)); }

        dim3_t dims() const { return dims_; }
        size_t size() const { return size_; }
        T * ptr() { return ptr_; }       
 
        T & operator[](int i) { return ptr_[i]; }
        T & operator()(int i, int j, int k){ 
            int idx = i * dims_.y * dims_.z + j * dims_.z + k;
            return ptr_[idx];
        }
        T & operator[](dim3_t i) { 
            int idx = i.x * dims_.y * dims_.z + i.y * dims_.z + i.z;
            return ptr_[idx];
        }
};

// pass by reference
template <typename T>
using Array_ref = Array<T> &;

#endif // ARRAY3__H
