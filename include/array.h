#include <algorithm>
#include <complex>
#include <cstdint>
#include <execution>
#include <memory>
#include <random>
#include <tuple>
#include <type_traits>

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

      public:
        Array() : dims_(0, 0, 0), size_(0), ptr_(nullptr) {}

        Array(size_t x, size_t y, size_t z) :
            dims_(x, y, z),
            size_(dims_.size()),
            ptr_(std::make_unique<T[]>(size_)) {}

        Array(dims_t d) :
            dims_(d), size_(dims_.size()), ptr_(std::make_unique<T[]>(size_)) {}

        Array(const Array<T> &) = delete;
        Array<T> &operator=(const Array<T> &) = delete;
        Array(Array<T> &&) noexcept = default;
        Array<T> &operator=(Array<T> &&) noexcept = default;

        [[nodiscard]] Array<T> clone() const {
            Array<T> rv(dims_);
            std::copy(this->begin(), this->end(), rv.begin());
            return std::move(rv);
        }

        [[nodiscard]] std::tuple<size_t, size_t, size_t> unravel_idx(
            size_t idx) {
            return dims_.unravel_idx(idx);
        }

        [[nodiscard]] size_t flatIdx(size_t i, size_t j, size_t k) const {
            return dims_.flat_idx(i, j, k);
        }

        [[nodiscard]] size_t flatIdx(int i0, int j0, int k0) const {
            auto i1 = static_cast<size_t>(i0);
            auto j1 = static_cast<size_t>(j0);
            auto k1 = static_cast<size_t>(k0);
            return dims_.flat_idx(i1, j1, k1);
        }

        void fill(T v) {
            std::fill(std::execution::par_unseq, this->begin(), this->end(), v);
        }

        void fill_random() {
            std::random_device rand_dev;
            std::mt19937 gen(rand_dev());
            std::uniform_real_distribution<> dist(0.0, 1.0);
            std::generate(this->begin(), this->end(),
                [&]() { return dist(gen); });
        }

        T *begin() { return ptr_.get(); }
        const T *begin() const { return ptr_.get(); }
        T *end() { return ptr_.get() + size_; }
        const T *end() const { return ptr_.get() + size_; }

        [[nodiscard]] dims_t dims() const { return dims_; }
        [[nodiscard]] size_t size() const { return size_; }
        [[nodiscard]] size_t nslices() const { return dims_.x(); }
        [[nodiscard]] size_t nrows() const { return dims_.y(); }
        [[nodiscard]] size_t ncols() const { return dims_.z(); }

        // indexing
        T &operator[](size_t i) { return ptr_[i]; }
        const T &operator[](size_t i) const { return ptr_[i]; }

    #if (__cplusplus == 202302L)
        T &operator[](size_t i, size_t j, size_t k) {
            return ptr_[flatIdx(i, j, k)];
        }
        T operator[](size_t i, size_t j, size_t k) const {
            return ptr_[flatIdx(i, j, k)];
        }
    #else
        T &operator[](dims_t i) { return ptr_[flatIdx(i.x(), i.y(), i.z())]; }
        const T &operator[](dims_t i) const {
            return ptr_[flatIdx(i.x(), i.y(), i.z())];
        }
    #endif

        // get slices
        T *slice(size_t i) { return ptr_.get() + (i * dims_.y() * dims_.z()); }
        const T *slice(size_t i) const {
            return ptr_.get() + (i * dims_.y() * dims_.z());
        }

        // multiplication operators
        Array<T> &operator*=(T v) {
            std::transform(std::execution::par_unseq, this->begin(),
                this->end(), this->begin(), [v](T x) { return x * v; });
            return *this;
        }
        Array<T> operator*(T v) {
            auto tmp = this->clone();
            tmp *= v;
            return tmp;
        }

        // division
        Array<T> &operator/=(T v) {
            if (std::abs(v) == 0) {
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
        Array<T> operator/(T scalar) const {
            auto rv = this->clone();
            rv /= scalar;
            return rv;
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
        Array<T> operator+(const Array<T> &rhs) const {
            auto tmp = this->clone();
            tmp += rhs;
            return tmp;
        }

        Array<T> &operator+=(const Array<T> &rhs) const {
            auto tmp = this->clone();
            tmp += rhs;
            return tmp;
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

        Array<T> operator-(const Array<T> &rhs) const {
            auto rv = this->clone();
            rv -= rhs;
            return rv;
        }

        static Array<T> zeros_like(const Array<T> &a) {
            Array<T> rv(a.dims());
            rv.fill(static_cast<T>(0));
            return rv;
        }
    };
} // namespace tomocam
#endif // ARRAY3__H
