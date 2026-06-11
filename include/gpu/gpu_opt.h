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

#ifndef TOMOCAM_GPU_OPT_H
#define TOMOCAM_GPU_OPT_H

#include <array>
#include <functional>

#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"

namespace tomocam::gpu::opt {

    template <typename T>
    class VecArray {
      private:
        std::array<DeviceArray<T>, 3> data_;

      public:
        VecArray() = default;
        VecArray(DeviceArray<T> a, DeviceArray<T> b, DeviceArray<T> c)
            : data_{std::move(a), std::move(b), std::move(c)} {}
        explicit VecArray(std::array<DeviceArray<T>, 3> &&arr)
            : data_(std::move(arr)) {}

        DeviceArray<T> &operator[](size_t i) { return data_[i]; }
        const DeviceArray<T> &operator[](size_t i) const { return data_[i]; }

        const std::array<DeviceArray<T>, 3> &data() const { return data_; }
        VecArray clone() const {
            return {data_[0].clone(), data_[1].clone(), data_[2].clone()};
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
            return {data_[0] + b[0], data_[1] + b[1], data_[2] + b[2]};
        }

        VecArray &operator+=(const VecArray &b) {
            data_[0] += b[0];
            data_[1] += b[1];
            data_[2] += b[2];
            return *this;
        }

        VecArray operator-(const VecArray &b) const {
            return {data_[0] - b[0], data_[1] - b[1], data_[2] - b[2]};
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
    using gpuFunction = std::function<VecArray<T>(const VecArray<T> &)>;

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

    /**
     * @brief GPU Conjugate Gradient solver for the preconditioned linear system
     *        A x = y, where A is a self-adjoint positive-definite operator.
     *
     * @param A        GPU operator (Compostion of RadonAdjoint and Radon)
     * @param y        Right-hand side: GPU array
     * @param x0       Initial guess
     * @param max_iter Maximum number of CG iterations
     * @param tol      Convergence tolerance (residual norm)
     * @return         Approximate solution x on the GPU
     */
    template <typename T>
    VecArray<T> cgsolver(const gpuFunction<T> &A, const VecArray<T> &yT,
                         const VecArray<T> &x0T, size_t max_iter, T tol, T xtol);

    /**
     * @brief GPU Split Bregman solver for the sparse angle laminography problem:
     *       min_x 0.5 ||A x - y||_2^2 + lambda * ||D x||_1
     *
     * @param A         GPU operator (Compostion of RadonAdjoint and Radon)
     * @param y         Backprojection of the measured data: GPU array
     * @param x0        Initial guess
     * @param lambda    Regularization parameter for the L1 term
     * @param mu        Bregman parameter (penalty parameter)
     * @param outer_max  Maximum number of outer iterations (Bregman updates)
     * @param inner_max  Maximum number of inner iterations (CG solver)
     * @param tol       Convergence tolerance for the CG solver (residual norm)
     * @param xtol      Convergence tolerance for the outer iterations (solution
     * change)
     * @return          Approximate solution x on the GPU
     */
    template <typename T>
    VecArray<T> split_bregman(const gpuFunction<T> &A, const VecArray<T> &y,
                              const VecArray<T> &x0, T lambda, T mu,
                              size_t outer_max, size_t inner_max, T tol, T xtol);

} // namespace tomocam::gpu::opt

#endif // TOMOCAM_GPU_OPT_H
