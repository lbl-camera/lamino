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
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "logger.h"

namespace tomocam::opt {

    template <typename T>
    using VecArray = std::array<Array<T>, 3>;

    template <typename T>
    VecArray<T> createVecArray() {
        return {Array<T>(), Array<T>(), Array<T>()};
    }

    template <typename T>
    VecArray<T> createVecArray(dims_t dims) {
        return {Array<T>(dims), Array<T>(dims), Array<T>(dims)};
    }

    template <typename T>
    VecArray<T> clone(const VecArray<T> &vec) {
        return {vec[0].clone(), vec[1].clone(), vec[2].clone()};
    }

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
    T dot(const VecArray<T> &a, const VecArray<T> &b) {
        T result = 0;
        for (size_t i = 0; i < 3; ++i) { result += array::dot(a[i], b[i]); }
        return result;
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
     * @param mask Mask for the support region
     * @param use_logfile Flag to enable logging
     * @return Reconstructed solution vector
     */
    template <typename T>
    VecArray<T> split_bregman(const Function<T> &A, const VecArray<T> &yT,
                              const VecArray<T> &x0, T lambda, T mu,
                              size_t outer_max, size_t inner_max, T tol, T xtol,
                              const Array<T> &mask, Logger *logger = nullptr);

    /**
     * @brief Constant-magnitude reconstruction with L2 (Tikhonov) smoothing,
     * solved by ADMM / scaled Bregman splitting:
     *   min_{x, M0} 0.5||Ax - y||^2 + lambda_l2/2 ||grad x||^2
     *   s.t. |x(r)| = M0 for r in the mask
     * with splitting variable z = x and dual b. The magnitude M0 is uniform but
     * unknown: it is initialised to the mean voxel magnitude of x0 and
     * re-estimated in every z-update as the mean magnitude of (x + b). The
     * constraint is only enforced inside the support mask; the remaining
     * voxels (padding) are left free.
     * @param A Function representing R^T R
     * @param yT Backprojected data R^T y
     * @param x0 Initial guess (also used to initialise M0)
     * @param lambda_l2 Weight of the L2 smoothness term
     * @param mu Penalty weight of the constraint splitting
     * @param outer_max Maximum number of outer iterations
     * @param inner_max Maximum number of CG iterations per outer iteration
     * @param tol Tolerance for the inner CG
     * @param xtol Tolerance on relative change of x and on ||x - z||/||x||
     * @param mask Binary support (same dims as x0, > 0.5 inside) where
     * |x| = M0 is enforced; pass all ones to constrain every voxel
     * @param zero_outside Restrict x and z to the mask: both are 0 outside it
     * (the CG solves the restricted system)
     * @return Solution z (|z| = M0 inside the mask)
     */
    template <typename T>
    VecArray<T> constmag_l2(const Function<T> &A, const VecArray<T> &yT,
                            const VecArray<T> &x0, T lambda_l2, T mu,
                            size_t outer_max, size_t inner_max, T tol, T xtol,
                            const Array<T> &mask, Logger *logger = nullptr,
                            bool zero_outside = false);

    /**
     * @brief Switchable regularization terms of regularized_admm. Each term is
     * independent: const_mag enables |x| = M0, lambda_tv > 0 enables TV and
     * lambda_l2 > 0 enables the L2 smoothness term.
     */
    template <typename T>
    struct RegWeights {
        bool const_mag = false;
        T lambda_tv = 0;
        T lambda_l2 = 0;

        [[nodiscard]] bool tv() const { return lambda_tv > T(0); }
        [[nodiscard]] bool l2() const { return lambda_l2 > T(0); }
    };

    /**
     * @brief Options of the magnitude constraint of regularized_admm. The
     * defaults reproduce the original behaviour (rho = mu, Ms re-estimated in
     * every iteration, no continuation, padding left free).
     */
    template <typename T>
    struct MagOptions {
        T rho = 0;                 // penalty of the z = x split; 0 means use mu
        T Ms = 0;                  // fixed magnitude if !estimate_Ms (0: mean |x0|)
        bool estimate_Ms = true;   // Ms <- mean_mask |x + b_z| before projecting
        T rho_growth = 1;          // rho <- min(rho * rho_growth, rho_max) per iter
        T rho_max = 0;             // 0 means unbounded
        bool zero_outside = false; // z = 0 outside the mask instead of leaving free
    };

    /**
     * @brief Per-iteration diagnostics of regularized_admm. Residuals are
     * divided by ||x||; terms that are switched
     * off are zero. data_shifted = 0.5<x,Ax> - <x,A^Ty>, i.e. the data misfit
     * without the constant 0.5||y||^2.
     */
    template <typename T>
    struct AdmmIter {
        size_t iter = 0;
        T primal_tv = 0, dual_tv = 0, primal_mag = 0, dual_mag = 0;
        T data_shifted = 0, tv = 0, l2 = 0;
        T rho = 0, Ms = 0;
    };

    template <typename T>
    struct AdmmHistory {
        std::vector<AdmmIter<T>> iters;
        VecArray<T> m;   // unprojected ADMM variable (|m|/Ms shows constraint load)
        VecArray<T> b_z; // final magnitude multiplier (scaled form)
        T rho = 0;       // final rho
        T Ms = 0;        // final Ms
    };

    /**
     * @brief Combined constant-magnitude / TV / L2 reconstruction, solved by
     * ADMM / scaled Bregman splitting:
     *   min_{x, M0} 0.5||Ax - y||^2 + lambda_tv TV(x) + lambda_l2/2 ||grad x||^2
     *   s.t. |x(r)| = M0 for r in the mask (only if const_mag)
     * with splitting variables z = x (magnitude, penalty rho) and d = grad x
     * (TV, penalty mu). Terms that are switched off cost no memory and no work;
     * TV only reproduces split_bregman and const_mag + L2 reproduces
     * constmag_l2. The sphere constraint is non-convex: only a stationary point
     * is reached, so warm-start from the unconstrained solution and use rho
     * continuation (mag.rho_growth > 1).
     * @param w Which terms are active and their weights
     * @param mu Penalty weight of the TV splitting (and of the magnitude
     * splitting unless mag.rho is set)
     * @param mask Binary support (same dims as x0, > 0.5 inside) where
     * |x| = M0 is enforced
     * @param mag Options of the magnitude constraint
     * @param history If not null, receives per-iteration diagnostics (this costs
     * one extra application of A per iteration)
     * @return z (|z| = M0 inside the mask) if const_mag, otherwise x
     */
    template <typename T>
    VecArray<T> regularized_admm(const Function<T> &A, const VecArray<T> &yT,
                                 const VecArray<T> &x0, const RegWeights<T> &w, T mu,
                                 size_t outer_max, size_t inner_max, T tol, T xtol,
                                 const Array<T> &mask, Logger *logger = nullptr,
                                 const MagOptions<T> &mag = MagOptions<T>(),
                                 AdmmHistory<T> *history = nullptr);

    /**
     * @brief Project a vector field onto the sphere |z| = M0 voxel-wise.
     * Zero vectors are left as zero.
     */
    template <typename T>
    VecArray<T> project_magnitude(const VecArray<T> &v, T M0);

    /**
     * @brief Voxel-wise projection onto |z| = Ms. Where |v| <= eps the
     * projection is not unique and the previous value z_prev is kept.
     */
    template <typename T>
    VecArray<T> project_sphere(const VecArray<T> &v, T Ms, const VecArray<T> &z_prev,
                               T eps);

    /**
     * @brief Mean voxel magnitude of a vector field.
     */
    template <typename T>
    T mean_magnitude(const VecArray<T> &v);

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
     * @param xtol Tolerance for convergence based on solution change (step norm)
     * @param mask Mask that defines support region
     * @param lambda Demagnetization constraint weight
     * @return Approximate solution vector
     */
    template <typename T>
    VecArray<T> cgsolver(const Function<T> &A, const VecArray<T> &b,
                         const VecArray<T> &x0, size_t max_iters, T tol, T xtol,
                         const Array<T> &mask, T lambda = 0,
                         Logger *logger = nullptr);
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
    void qggmrf(Array<T> &g, const Array<T> &x, T sigma, T p);

} // namespace tomocam::opt

#endif // TOMOCAM_OPT_H
