
#ifndef ARRAY3__H
#define ARRAY3__H

#include <string.h>

#include "dtypes.h"

template <typename T>
class Array {
 private:
    dim2_t dims_;
    size_t size_;
    T *ptr_;

 public:
    Array(int x, int y): dims_{x, y}, size_(x * y), ptr_(new T[x * y]) {}
    ~Array() { if (ptr_) delete[]ptr_; }

    Array(const Array & rhs) {
        dims_ = rhs.dims_;
        size_ = rhs.size_;
        ptr_ = new T [size_];
        memcpy(ptr_, rhs.ptr_, size_ * sizeof(T));
    }

    Array & operator=(const Array & rhs) {
        dims_ = rhs.dims_;
        size_ = rhs.size_;
        ptr_ = new T [size_];
        memcpy(ptr_, rhs.ptr_, size_ * sizeof(T));
        return *this;
    }
     
    void clear() {
        memset(ptr_, 0, size_ * sizeof(T));
    }

    dim2_t dims() const {
        return dims_;
    }
    size_t size() const {
        return size_;
    }
    T *ptr() {
        return ptr_;
    }

    T & operator[](int i) { return ptr_[i]; }

    T & operator()(int i, int j) {
        int idx = i * dims_.y + j;
        return ptr_[idx];
    }

    T & operator[](dim2_t i) {
        int idx = i.x * dims_.y + i.y;
        return ptr_[idx];
    }

	Array<float> real() {
		Array<float> rv(dims_.x, dims_.y);
		for (int i = 0; i < size_; i++)
			rv[i] = ptr_[i].real();
		return rv;
	}
};

#endif // ARRAY3__H
