#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <future>
#include <memory>
#include <string.h>
#include <tuple>

#include "dtypes.h"

#ifndef ARRAY3__H
#define ARRAY3__H

namespace tomocam {
    template <typename T> class Array {
      private:
        dims_t dims_;
        uint64_t size_;
        std::unique_ptr<T[]> ptr_;

      public:
        Array() : dims_(0, 0, 0), size_(0), ptr_(nullptr) {}

        Array(uint64_t x, uint64_t y, uint64_t z)
            : dims_(x, y, z), size_(dims_.size()), ptr_(std::make_unique<T[]>(size_)) {}

        Array(dims_t d) : dims_(d), size_(dims_.size()), ptr_(std::make_unique<T[]>(size_)) {}

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

        std::tuple<uint64_t, uint64_t, uint64_t> unravel_idx(uint64_t idx) {
            return dims_.unravel_idx(idx);
        }

        uint64_t flatIdx(uint64_t i, uint64_t j, uint64_t k) const {
            return dims_.flat_idx(i, j, k);
        }

        uint64_t flatIdx(int i0, int j0, int k0) const {
            uint64_t i1 = static_cast<uint64_t>(i1);
            uint64_t j1 = static_cast<uint64_t>(j1);
            uint64_t k1 = static_cast<uint64_t>(k1);
            return dims_.flat_idx(i1, j1, k1);
        }

        void fill(T v) { std::fill(std::execution::par_unseq, this->begin(), this->end(), v); }

        T *begin() { return ptr_.get(); }
        const T *begin() const { return ptr_.get(); }
        T *end() { return ptr_.get() + size_; }
        const T *end() const { return ptr_.get() + size_; }

        dims_t dims() const { return dims_; }
        uint64_t size() const { return size_; }
        uint64_t nslices() const { return dims_.x(); }
        uint64_t nrows() const { return dims_.y(); }
        uint64_t ncols() const { return dims_.z(); }

        // indexing
        T &operator[](uint64_t i) { return ptr_[i]; }
        T operator[](uint64_t i) const { return ptr_[i]; }

        T &operator[](uint64_t i, uint64_t j, uint64_t k) { return ptr_[flatIdx(i, j, k)]; }
        T operator[](uint64_t i, uint64_t j, uint64_t k) const { return ptr_[flatIdx(i, j, k)]; }

        T &operator[](dims_t i) { return ptr_[flatIdx(i.x(), i.y(), i.z())]; }
        const T &operator[](dims_t i) const { return ptr_[flatIdx(i.x(), i.y(), i.z())]; }

        // multiplication operators
        Array<T> &operator*=(T v) {
            std::transform(std::execution::par_unseq, this->begin(), this->end(), ptr_.get(),
                           [v](T x) { return x * v; });
            return *this;
        }

        // division
        Array<T> &operator/=(T v) {
            if (v == 0) {
                throw std::runtime_error("Division by zero in Array<T>::operator/=");
            }
            std::transform(std::execution::par_unseq, this->begin(), this->end(), ptr_.get(),
                           [v](T x) { return x / v; });
            return *this;
        }
        Array<T> &operator/=(const Array<T> &v) {
            std::transform(std::execution::par_unseq, this->begin(), this->end(), v.ptr_.get(),
                           ptr_.get(), [](T x, T y) {
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
            std::transform(std::execution::par_unseq, this->begin(), this->end(), ptr_.get(),
                           [v](T x) { return x + v; });
            return *this;
        }
        Array<T> &operator+=(const Array<T> &v) {
            std::transform(std::execution::par_unseq, this->begin(), this->end(), v.ptr_.get(),
                           ptr_.get(), std::plus<T>());
            return *this;
        }

        // subtraction
        Array<T> &operator-=(T v) {
            std::transform(std::execution::par_unseq, this->begin(), this->end(), ptr_.get(),
                           [v](T x) { return x - v; });
            return *this;
        }
        Array<T> &operator-=(const Array<T> &v) {
            std::transform(std::execution::par_unseq, this->begin(), this->end(), v.ptr_.get(),
                           ptr_.get(), std::minus<T>());
            return *this;
        }
    };

    // multiplications
    template <typename T> Array<T> &operator*(Array<T> a, const Array<T> &b) { return a *= b; }
    template <typename T> Array<T> &operator*(T a, Array<T> b) { return b *= a; }
    template <typename T> Array<T> &operator*(Array<T> b, T a) { return b *= a; }
    //  addtitions
    template <typename T> Array<T> &operator+(Array<T> a, const Array<T> &b) { return a += b; }
    template <typename T> Array<T> &operator+(T a, Array<T> b) { return b += a; }
    template <typename T> Array<T> &operator+(Array<T> b, T a) { return b += a; }
    // subtractions
    template <typename T> Array<T> &operator-(Array<T> a, const Array<T> &b) { return a -= b; }

    template <typename T> Array<T> operator-(double a, Array<T> b) {
        Array<T> rv(b.dims());
        std::transform(std::execution::par_unseq, b.begin(), b.end(), rv.begin(),
                       [a](T x) { return a - x; });
        return rv;
    }

    // divisions
    template <typename T> Array<T> &operator/(Array<T> a, const Array<T> &b) { return a /= b; }
    template <typename T> Array<T> &operator/(Array<T> a, T b) { return a /= b; }
    template <typename T> Array<T> operator/(double a, Array<T> b) {
        Array<T> rv(b.dims());
        std::transform(std::execution::par_unseq, b.begin(), b.end(), rv.begin(),
                       [a](T x) { return a / x; });
        return rv;
    }

} // namespace tomocam
#endif // ARRAY3__H
