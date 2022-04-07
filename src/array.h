
#ifndef ARRAY3__H
#define ARRAY3__H

#include <string.h>
#include <fstream>

#include "dtypes.h"

template <typename T>
class Array {
 private:
    dim2_t dims_;
    size_t size_;
    T *ptr_;

 public:
    Array(int x, int y): dims_{x, y}, size_(x * y), ptr_(new T[x * y]) {}
    Array(int x, int y, T v): dims_{x, y}, size_(x * y), ptr_(new T[x * y]) {
        for (int i = 0; i < x*y; i++)
            ptr_[i] = v;
    }
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

    dim2_t dims() const { return dims_; }
    size_t size() const { return size_; }
    T *ptr() { return ptr_; }
    const T *ptr() const { return ptr_; }

    // indexing
    T & operator[](int i) { return ptr_[i]; }
    const T & operator[](int i) const { return ptr_[i]; }

    T & operator()(int i, int j) { return ptr_[i * dims_.y + j]; }
    const T & operator()(int i, int j) const { return ptr_[i * dims_.y + j]; }

    T & operator[](dim2_t i) { return ptr_[i.x * dims_.y + i.y]; }
    const T & operator[](dim2_t i) const { return ptr_[i.x * dims_.y + i.y]; }

    // convert to complex
    Array<complex_t> cmplx() {
        Array<complex_t> v(dims_.x, dims_.y);
        #pragma omp parllel for
        for (int i = 0; i < size_; i++)
            v[i] = static_cast<complex_t>(ptr_[i]);
        return v;
    }

    // multiplication operators
    Array<T>& operator*=(T v) {
        for (int i = 0; i < size_; i++) 
            ptr_[i] *= v;
        return *this;
    }
    Array<T>& operator*=(const Array<T> &v) {
        for (int i = 0; i < size_; i++) 
            ptr_[i] *= v.ptr_[i];
        return *this;
    }

    // division
    Array<T>& operator/=(T v) {
        for (int i = 0; i < size_; i++) 
            ptr_[i] /= v;
        return *this;
    }
    Array<T>& operator/=(const Array<T> &v) {
        for (int i = 0; i < size_; i++) 
            ptr_[i] /= v.ptr_[i];
        return *this;
    }

    // addition
    Array<T>& operator+=(T v) {
        for (int i = 0; i < size_; i++) 
            ptr_[i] += v;
        return *this;
    }
    Array<T>& operator+=(const Array<T> &v) {
        for (int i = 0; i < size_; i++) 
            ptr_[i] += v.ptr_[i];
        return *this;
    }

    // subtraction
    Array<T>& operator-=(T v) {
        for (int i = 0; i < size_; i++) 
            ptr_[i] -= v;
        return *this;
    }
    Array<T>& operator-=(const Array<T> &v) {
        for (int i = 0; i < size_; i++) 
            ptr_[i] -= v.ptr_[i];
        return *this;
    }

    Array<T> sqrt() const {
        Array<T> rv(dims_.x, dims_.y);
        for (int i = 0; i < size_; i++)
            rv.ptr_[i] = std::sqrt(ptr_[i]);
        return rv;
    } 

    // get real part
	Array<float> real() {
		Array<float> rv(dims_.x, dims_.y);
		for (int i = 0; i < size_; i++)
			rv[i] = ptr_[i].real();
		return rv;
	}

    // norm2
    float norm() const {
        float v = 0;
        for (int i = 0; i < size_; i++)
            v += std::pow(std::abs(ptr_[i]), 2);
        return std::sqrt(v);
    }

    // max
    T max() const {
        T maxv = 1 << 30;
        for (int i = 0; i < size_; i++)
            if (ptr_[i] > maxv) 
                maxv = ptr_[i];
        return maxv;
    } 
        
    // read binary file
    void fromfile(const char *filename) const {
        std::ifstream in(filename, std::ios::binary | std::ios::in);
        in.read((char *) ptr_, size_ * sizeof(T));
        in.close();
    }
 
    // write to binary file
    void tofile(const char *filename) const {
        std::ofstream out(filename, std::ios::binary | std::ios::in);
        out.write((char *) ptr_, size_ * sizeof(T));
        out.close();
    }

    // transpose
    void transpose() { 
        for (int i = 1; i < dims_.x; i++)
            for (int j = 0; j < i; j++) {
                auto temp  = ptr_[i * dims_.y + j];
                ptr_[i * dims_.y + j] = ptr_[i + j * dims_.x]; 
                ptr_[i + dims_.x * j] = temp;
            }
    }
};

// multiplications
template <typename T>
Array<T>& operator*(Array<T> a, const Array<T> &b) {
    return a *= b;
}
template <typename T>
Array<T>& operator*(T a, Array<T> b) {
    return b *= a;
}
template <typename T>
Array<T>& operator*(Array<T> b, T a) {
    return b *= a;
}
//  addtitions
template <typename T>
Array<T>& operator+(Array<T> a, const Array<T> &b) {
    return a += b; 
}
template <typename T>
Array<T>& operator+(T a, Array<T> b) {
    return b += a;
}
template <typename T>
Array<T>& operator+(Array<T> b, T a) {
    return b += a;
}
// subtractions
template <typename T>
Array<T>& operator-(Array<T> a, const Array<T> &b) {
    return a -= b;
}
template <typename T>
Array<T>& operator-(float a, Array<T> b) {
    return b -= a;
}

// divisions
template <typename T>
Array<T>& operator/(Array<T> a, const Array<T> &b) {
    return a /= b;
}
template <typename T>
Array<T>& operator/(Array<T> a, T b) {
    return a /= b;
}
template <typename T>
Array<T>& operator/(float a, Array<T> b) {
    Array<T> rv(b.dims().x, b.dims().y);
    for (int i = 0; i < b.size(); i++)
        rv[i] = a / b[i];
    return rv;
}

#endif // ARRAY3__H
