/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 * National Laboratory (subject to receipt of any required approvals from the
 * U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 * IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 * the U.S. Government has been granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 * to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#ifndef VEC_ARRAY_H
#define VEC_ARRAY_H

#include <array>

#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"

namespace tomocam::gpu {

    template <typename T>
    class VecArray {
      private:
        std::array<DeviceArray<T>, 3> data_;

      public:
        VecArray() = default;
        explicit VecArray(DeviceArray<T> &&x, DeviceArray<T> &&y, DeviceArray<T> &&z)
            : data_{std::move(x), std::move(y), std::move(z)} {}
        explicit VecArray(std::array<DeviceArray<T>, 3> &&data)
            : data_(std::move(data)) {}

        DeviceArray<T> &operator[](size_t i) { return data_[i]; }
        const DeviceArray<T> &operator[](size_t i) const { return data_[i]; }

        const std::array<DeviceArray<T>, 3> &data() const { return data_; }
        VecArray clone() const {
            return VecArray{data_[0].clone(), data_[1].clone(), data_[2].clone()};
        }
        static VecArray zeros(dims_t dims) {
            return {DeviceArray<T>::zeros(dims), DeviceArray<T>::zeros(dims),
                    DeviceArray<T>::zeros(dims)};
        }

        T dot(const VecArray &b) const {
            return array::dot(data_[0], b[0]) + array::dot(data_[1], b[1]) +
                   array::dot(data_[2], b[2]);
        }

        T norm2() const { return sqrt(dot(*this)); }

        VecArray operator+(const VecArray &b) const {
            return VecArray{data_[0] + b[0], data_[1] + b[1], data_[2] + b[2]};
        }

        VecArray &operator+=(const VecArray &b) {
            data_[0] += b[0];
            data_[1] += b[1];
            data_[2] += b[2];
            return *this;
        }

        VecArray operator-(const VecArray &b) const {
            return VecArray{data_[0] - b[0], data_[1] - b[1], data_[2] - b[2]};
        }

        VecArray &operator-=(const VecArray &b) {
            data_[0] -= b[0];
            data_[1] -= b[1];
            data_[2] -= b[2];
            return *this;
        }

        VecArray operator*(T s) const {
            return {data_[0] * s, data_[1] * s, data_[2] * s};
        }

        VecArray &operator*=(T s) {
            data_[0] *= s;
            data_[1] *= s;
            data_[2] *= s;
            return *this;
        }
    };

    template <typename T>
    void vec_xpay(VecArray<T> &x, const VecArray<T> &y, T alpha) {
        array::xpay(x[0], y[0], alpha);
        array::xpay(x[1], y[1], alpha);
        array::xpay(x[2], y[2], alpha);
    }

    template <typename T>
    void vec_axpy(VecArray<T> &x, T alpha, const VecArray<T> &y) {
        array::axpy(x[0], alpha, y[0]);
        array::axpy(x[1], alpha, y[1]);
        array::axpy(x[2], alpha, y[2]);
    }
} // namespace tomocam::gpu

#endif // VEC_ARRAY_H
