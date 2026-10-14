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
#include <vector>

#include "dtypes.h"
#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/vec_array.h"
#include "logger.h"

namespace tomocam::gpu::opt {

    template <typename T>
    using gpuFunction = std::function<VecArray<T>(const VecArray<T> &)>;

    /**
     * @brief GPU Conjugate Gradient solver for the preconditioned linear system
     *        A x = y, where A is a self-adjoint positive-definite operator.
     *
     * @param A        GPU operator (Compostion of RadonAdjoint and Radon)
     * @param y        Right-hand side: GPU array
     * @param x0       Initial guess
     * @param max_iter Maximum number of CG iterations
     * @param tol      Convergence tolerance (residual norm)
     * @param support  Binary support: x is kept 0 where it is 0 (CG on the
     *                 restricted system); an empty array means no constraint
     * @return         Approximate solution x on the GPU
     */
    template <typename T>
    VecArray<T> cgsolver(const gpuFunction<T> &A, const VecArray<T> &yT,
                         const VecArray<T> &x0T, size_t max_iter, T tol, T xtol,
                         const DeviceArray<T> &support, T lambda = 0,
                         Logger *logger = nullptr);

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
                              size_t outer_max, size_t inner_max, T tol, T xtol,
                              const DeviceArray<T> &support,
                              Logger *logger = nullptr);

    /**
     * @brief GPU constant-magnitude reconstruction with L2 (Tikhonov) smoothing,
     * solved by ADMM / scaled Bregman splitting (mirrors
     * tomocam::opt::constmag_l2):
     *   min_{x, M0} 0.5||Ax - y||^2 + lambda_l2/2 ||grad x||^2
     *   s.t. |x(r)| = M0 for r in the mask
     * The magnitude M0 is uniform but unknown: it is initialised to the mean
     * voxel magnitude of x0 and re-estimated in every z-update. The constraint
     * is only enforced inside the support mask; the remaining voxels (padding)
     * are left free.
     *
     * @param A         GPU operator (Compostion of RadonAdjoint and Radon)
     * @param yT        Backprojection of the measured data: GPU array
     * @param x0        Initial guess (also used to initialise M0)
     * @param lambda_l2 Weight of the L2 smoothness term
     * @param mu        Penalty weight of the constraint splitting
     * @param outer_max Maximum number of outer iterations
     * @param inner_max Maximum number of CG iterations per outer iteration
     * @param tol       Tolerance for the inner CG
     * @param xtol      Tolerance on relative change of x and on ||x - z||/||x||
     * @param mask      Binary support (same dims as x0, > 0.5 inside) where
     * |x| = M0 is enforced; all ones constrains every voxel
     * @return          Solution z on the GPU (|z| = M0 inside the mask)
     */
    template <typename T>
    VecArray<T> constmag_l2(const gpuFunction<T> &A, const VecArray<T> &yT,
                            const VecArray<T> &x0, T lambda_l2, T mu,
                            size_t outer_max, size_t inner_max, T tol, T xtol,
                            const DeviceArray<T> &mask, Logger *logger = nullptr,
                            bool zero_outside = false);

    /// Switchable terms of regularized_admm (see tomocam::opt::RegWeights)
    template <typename T>
    struct RegWeights {
        bool const_mag = false;
        T lambda_tv = 0; ///< > 0 enables TV
        T lambda_l2 = 0; ///< > 0 enables the L2 smoothness term

        bool tv() const { return lambda_tv > T(0); }
        bool l2() const { return lambda_l2 > T(0); }
    };

    /// Options of the magnitude constraint (see tomocam::opt::MagOptions)
    template <typename T>
    struct MagOptions {
        T rho = 0;                 ///< penalty of the z = x split; 0 means use mu
        T Ms = 0;                  ///< fixed magnitude if !estimate_Ms
        bool estimate_Ms = true;   ///< Ms <- mean_mask |x + b_z| before projecting
        T rho_growth = 1;          ///< rho <- min(rho * rho_growth, rho_max)
        T rho_max = 0;             ///< 0 means unbounded
        bool zero_outside = false; ///< z = 0 outside the mask instead of free
    };

    /// Per-iteration diagnostics. Unlike the CPU solver only the magnitude
    /// residuals, rho and Ms are recorded (no TV residuals or objective terms).
    template <typename T>
    struct AdmmIter {
        size_t iter = 0;
        T primal_mag = 0, dual_mag = 0;
        T rho = 0, Ms = 0;
    };

    template <typename T>
    struct AdmmHistory {
        std::vector<AdmmIter<T>> iters;
        VecArray<T> m;   ///< last unprojected variable (const_mag only)
        VecArray<T> b_z; ///< final magnitude multiplier (const_mag only)
        T rho = 0;
        T Ms = 0;
    };

    /**
     * @brief GPU combined constant-magnitude / TV / L2 reconstruction (mirrors
     * tomocam::opt::regularized_admm):
     *   min_{x, M0} 0.5||Ax - y||^2 + lambda_tv TV(x) + lambda_l2/2 ||grad x||^2
     *   s.t. |x(r)| = M0 for r in the mask (only if const_mag)
     * Terms that are switched off allocate nothing and cost no work.
     * @return z (|z| = M0 inside the mask) if const_mag, otherwise x
     */
    template <typename T>
    VecArray<T> regularized_admm(const gpuFunction<T> &A, const VecArray<T> &yT,
                                 const VecArray<T> &x0, const RegWeights<T> &w, T mu,
                                 size_t outer_max, size_t inner_max, T tol, T xtol,
                                 const DeviceArray<T> &mask,
                                 Logger *logger = nullptr,
                                 const MagOptions<T> &mag = MagOptions<T>(),
                                 AdmmHistory<T> *history = nullptr);

    /**
     * @brief Voxel-wise projection onto |z| = Ms; where |v| <= eps the previous
     * value z_prev is kept.
     */
    template <typename T>
    VecArray<T> project_sphere(const VecArray<T> &v, T Ms, const VecArray<T> &z_prev,
                               T eps);

    /**
     * @brief Project a vector field onto the sphere |z| = M0 voxel-wise.
     * Zero vectors are left as zero.
     */
    template <typename T>
    VecArray<T> project_magnitude(const VecArray<T> &v, T M0);

    /**
     * @brief Mean voxel magnitude of a vector field.
     */
    template <typename T>
    T mean_magnitude(const VecArray<T> &v);

} // namespace tomocam::gpu::opt

#endif // TOMOCAM_GPU_OPT_H
