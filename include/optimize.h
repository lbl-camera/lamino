/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 *National Laboratory (subject to receipt of any required approvals from the
 *U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 *IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 *the U.S. Government has been granted for itself and others acting on its
 *behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 *to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#ifndef TOMOCAM_OPT_H
#define TOMOCAM_OPT_H

#include <cmath>
#include <iostream>

#include "array.h"
#include "array_ops.h"
#include "logger.h"

namespace tomocam::opt {

    template <typename T>
    using VecArray = std::array<Array<T>, 3>;

    // addtion operator for VecArray
    template <typename T>
    VecArray<T> operator+(const VecArray<T> &a, const VecArray<T> &b) {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    }
    // in-place addition operator for VecArray
    template <typename T>
    VecArray<T> &operator+=(VecArray<T> &a, const VecArray<T> &b) {
        a[0] += b[0];
        a[1] += b[1];
        a[2] += b[2];
        return a;
    }
    // subtraction operator for VecArray
    template <typename T>
    VecArray<T> operator-(const VecArray<T> &a, const VecArray<T> &b) {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    }
    // in-place subtraction operator for VecArray
    template <typename T>
    VecArray<T> &operator-=(VecArray<T> &a, const VecArray<T> &b) {
        a[0] -= b[0];
        a[1] -= b[1];
        a[2] -= b[2];
        return a;
    }
    // scalar multiplication operator for VecArray
    template <typename T>
    VecArray<T> operator*(const VecArray<T> &a, T scalar) {
        return {a[0] * scalar, a[1] * scalar, a[2] * scalar};
    }

    template <typename T>
    using Function = std::function<VecArray<T>(const VecArray<T> &)>;

    template <typename T>
    using Residual = std::function<T(const VecArray<T> &)>;

    /** implement Split Bregman method for L1-regularized optimization problems
     * @brief Split Bregman method for L1-regularized optimization problems
     * @param A Function representing the R^T R operator for datasets 1, 2, ... (i.e.
     * A(u) = [R1^T R1 u, R2^T R2 u, ...])
     * @param yT 3-component arrays representing the R1^T y1 + R2^T y2 + ...
     * @param x0 Initial guess for the solution (std::array of 3 components)
     * @param lambda Regularization parameter
     * @param mu Bregman parameter
     * @param outer_max Maximum number of outer iterations
     * @param inner_max Maximum number of inner iterations
     * @param tol Tolerance for convergence based on residual norm
     * @param xtol Tolerance for convergence based on solution change
     * @return Reconstructed solution vector
     */
    template <typename T>
    VecArray<T> split_bregman(const Function<T> &A, const VecArray<T> &yT,
                              const VecArray<T> &x0, T lambda, T mu,
                              size_t outer_max, size_t inner_max, T tol, T xtol);

    /**
     * @brief Conjugate Gradient method for solving the
     * linear system [A1, A2, ...]^T x = [b1, b2, ...]^T, where A1, A2 are
     * self-adjoint operators.
     *
     * @param A(u) Function represeting the combined (A1 + A2 + ... + μ∇^T∇)u
     * @param b 3-component array (R^T y1 + * R^T y2 + ... + μ∇^T(d- b))
     * @param x0 Initial guess for the solution
     * @param max_iters Maximum number of iterations
     * @param tol Tolerance for convergence based on residual norm
     * @param lambda Demagnetization constraint weight
     * @return Approximate solution vector
     */
    template <typename T>
    VecArray<T> cgsolver(const Function<T> &A, const VecArray<T> &b,
                         const VecArray<T> &x0, size_t max_iters, T tol, T lambda);
    /**
     * @brief Nesterov's Optimal Gradient Method with Boyd's momentum term (vector
     * version)
     * @param grad Gradient function for vector reconstruction
     * @param loss Loss function for vector reconstruction
     * @param x Initial solution (array of 3 components)
     * @param max_iters Maximum number of iterations
     * @param lipschitz Lipschitz constant of the gradient
     * @param tol Tolerance for convergence based on loss change
     * @param xtol Tolerance for convergence based on solution change
     * @param max_inner_iters Maximum number of inner iterations for line search
     * @param logger Logger instance for output control
     * @return Optimized solution (array of 3 components)
     */
    template <typename T>
    std::array<Array<T>, 3> nagopt(
        const std::function<std::array<Array<T>, 3>(const std::array<Array<T>, 3> &)>
            &grad,
        const std::function<T(const std::array<Array<T>, 3> &)> &loss,
        std::array<Array<T>, 3> &x, size_t max_iters, T lipschitz, T tol, T xtol,
        size_t max_inner_iters = 20, Logger *logger = nullptr);

    /**
     * @brief Estimate the Lipschitz constant of a gradient function using the power
     * iteration method
     * @param function representing the system matrix
     * @param x0 reference for the size of the input
     * @param max_iters Maximum number of iterations (default: 20)
     * @param tol Tolerance for convergence (default: 1e-5)
     * @return Estimated Lipschitz constant
     */
    template <typename T>
    T lipschitz(const Function<T> &grad, const Array<T> &x0, size_t max_iters = 20,
                T tol = 1e-5);

    /**
     * @brief Compute the gradient of the q-generalized Gaussian Markov Random Field
     * (qGGMRF) prior
     * @param x Input array
     * @param dx Output gradient array (updated in place)
     * @param sigma Scale parameter
     * @param p Exponent parameter
     */
    template <typename T>
    void qggmrf(const Array<T> &x, Array<T> &dx, T sigma, T p);

} // namespace tomocam::opt

#endif // TOMOCAM_OPT_H
