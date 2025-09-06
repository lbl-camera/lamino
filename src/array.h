#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <memory>
#include <string.h>
#include <tuple>

#include "dtypes.h"

#ifndef ARRAY3__H
#define ARRAY3__H

namespace tomocam {
    template <typename T>
    class Array {
      private:
        dims_t dims_;
        size_t size_;
        std::unique_ptr<T[]> ptr_;

        size_t flatIdx(int i, int j, int k) {
            return (static_cast<size_t>(i) * dims_.x + static_cast<size_t>(j)) *
                       dims_.y +
                   static_cast<size_t>(k);
        }
        size_t mul(int i, int j, int k) {
            return (static_cast<size_t>(i) * static_cast<size_t>(j) *
                    static_cast<size_t>(k));
        }

      public:
        Array() : dims_{0, 0, 0}, size_(0), ptr_(nullptr) {}

        Array(int x, int y, int z) :
            dims_(x, y, z),
            size_(mul(x, y, z)),
            ptr_(std::make_unique<T[]>(mul(x, y, z))) {}

        Array(int x, int y, int z, T v) :
            dims_(x, y, z),
            size_(mul(x, y, z)),
            ptr_(std::make_unique<T[]>(mul(x, y, z))) {
            std::fill(std::execution::par_unseq, ptr_.get(), ptr_.get() + size_,
                v);
        }

        Array(dims_t d) :
            dims_(d),
            size_(mul(d.x, d.y, d.z)),
            ptr_(std::make_unique<T[]>(mul(d.x, d.y, d.z))) {}

        Array(const Array &rhs) {
            dims_ = rhs.dims_;
            size_ = rhs.size_;
            ptr_ = std::make_unique<T[]>(size_);
            std::copy(rhs.begin(), rhs.end(), ptr_.get());
        }

        Array &operator=(const Array &rhs) {
            if (this != &rhs) {
                dims_ = rhs.dims_;
                size_ = rhs.size_;
                ptr_.reset();
                ptr_ = std::make_unique<T[]>(size_);
                std::copy(rhs.begin(), rhs.end(), ptr_.get());
            }
            return *this;
        }

        Array &operator=(T v) {
            std::fill(std::execution::par_unseq, this->begin(), this->end(), v);
            return *this;
        }

        std::tuple<size_t, size_t, size_t> unravel_idx(size_t idx) {
            size_t i = idx / dims_.y / dims_.z;
            size_t j = (idx / dims_.z) % dims_.y;
            size_t k = idx % dims_.z;
            return std::make_tuple(i, j, k);
        }

        size_t flat_idx(size_t i, size_t j, size_t k) {
            return (i * dims_.y * dims_.z + j * dims_.z + k);
        }

        T *begin() { return ptr_.get(); }
        const T *begin() const { return ptr_.get(); }
        T *end() { return ptr_.get() + size_; }
        const T *end() const { return ptr_.get() + size_; }

        dims_t dims() const { return dims_; }
        size_t size() const { return size_; }
        int nslices() const { return dims_.x; }
        int nrows() const { return dims_.y; }
        int ncols() const { return dims_.z; }

        // indexing
        T &operator[](int i) { return ptr_[i]; }
        T operator[](int i) const { return ptr_[i]; }

        T &operator[](int i, int j, int k) { return ptr_[flatIdx(i, j, k)]; }
        T operator()(int i, int j, int k) const {
            return ptr_[flatIdx(i, j, k)];
        }

        T &operator[](dims_t i) { return ptr_[flatIdx(i.x, i.y, i.z)]; }
        T operator[](dims_t i) const { return ptr_[flatIdx(i.x, i.y, i.z)]; }

        // multiplication operators
        Array<T> &operator*=(T v) {
            std::transform(std::execution::par_unseq, this->begin(),
                this->end(), ptr_.get(), [v](T x) { return x * v; });
            return *this;
        }

        // division
        Array<T> &operator/=(T v) {
            if (v == 0) {
                throw std::runtime_error(
                    "Division by zero in Array<T>::operator/=");
            }
            std::transform(std::execution::par_unseq, this->begin(),
                this->end(), ptr_.get(), [v](T x) { return x / v; });
            return *this;
        }
        Array<T> &operator/=(const Array<T> &v) {
            std::transform(std::execution::par_unseq, this->begin(),
                this->end(), v.ptr_.get(), ptr_.get(), [](T x, T y) {
                    if (y == 0) {
                        throw std::runtime_error(
                            "Division by zero in Array<T>::operator/=");
                    }
                    return x / y;
                });
            return *this;
        }

        // addition
        Array<T> &operator+=(T v) {
            std::transform(std::execution::par_unseq, this->begin(),
                this->end(), ptr_.get(), [v](T x) { return x + v; });
            return *this;
        }
        Array<T> &operator+=(const Array<T> &v) {
            std::transform(std::execution::par_unseq, this->begin(),
                this->end(), v.ptr_.get(), ptr_.get(), std::plus<T>());
            return *this;
        }

        // subtraction
        Array<T> &operator-=(T v) {
            std::transform(std::execution::par_unseq, this->begin(),
                this->end(), ptr_.get(), [v](T x) { return x - v; });
            return *this;
        }
        Array<T> &operator-=(const Array<T> &v) {
            std::transform(std::execution::par_unseq, this->begin(),
                this->end(), v.ptr_.get(), ptr_.get(), std::minus<T>());
            return *this;
        }
    };

    // multiplications
    template <typename T>
    Array<T> &operator*(Array<T> a, const Array<T> &b) {
        return a *= b;
    }
    template <typename T>
    Array<T> &operator*(T a, Array<T> b) {
        return b *= a;
    }
    template <typename T>
    Array<T> &operator*(Array<T> b, T a) {
        return b *= a;
    }
    //  addtitions
    template <typename T>
    Array<T> &operator+(Array<T> a, const Array<T> &b) {
        return a += b;
    }
    template <typename T>
    Array<T> &operator+(T a, Array<T> b) {
        return b += a;
    }
    template <typename T>
    Array<T> &operator+(Array<T> b, T a) {
        return b += a;
    }
    // subtractions
    template <typename T>
    Array<T> &operator-(Array<T> a, const Array<T> &b) {
        return a -= b;
    }

    template <typename T>
    Array<T> operator-(double a, Array<T> b) {
        Array<T> rv(b.dims());
        std::transform(std::execution::par_unseq, b.begin(), b.end(),
            rv.begin(), [a](T x) { return a - x; });
        return rv;
    }

    // divisions
    template <typename T>
    Array<T> &operator/(Array<T> a, const Array<T> &b) {
        return a /= b;
    }
    template <typename T>
    Array<T> &operator/(Array<T> a, T b) {
        return a /= b;
    }
    template <typename T>
    Array<T> operator/(double a, Array<T> b) {
        Array<T> rv(b.dims());
        std::transform(std::execution::par_unseq, b.begin(), b.end(),
            rv.begin(), [a](T x) { return a / x; });
        return rv;
    }

} // namespace tomocam
#endif // ARRAY3__H
