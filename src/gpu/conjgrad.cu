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

#include <cstdio>
#include <format>
#include <limits>

#include <cuda_profiler_api.h>

#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/transform.h>

#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/gpu_opt.h"
#include "gpu/mem_check.h"
#include "gpu/utils.h"
#include "gpu/vec_array.h"
#include "mask.h"

namespace tomocam::gpu::opt {

    template <typename T>
    VecArray<T> cgsolver(const gpuFunction<T> &A, const VecArray<T> &y,
                         const VecArray<T> &x0, size_t max_iter, T tol, T xtol,
                         dims_t support_dims, T lambda, Logger *logger) {

        auto precond_apply = [](const VecArray<T> &r) {
            return r.clone(); // placeholder
        };

        // Build support mask and upload to GPU
        auto cpu_mask = mask_support<T>(x0[0].dims(), support_dims);
        DeviceArray<T> gpu_mask(cpu_mask);
        auto apply_support = [&gpu_mask](VecArray<T> &v) {
            for (size_t i = 0; i < 3; ++i) { v[i] *= gpu_mask; }
        };

        // Initialize solution and residual arrays
        VecArray<T> x = x0.clone();

        // r = y - A(x)
        auto r = y - A(x);
        // z = M^{-1} r,  p = z,  rs_old = z^T r
        VecArray<T> z = precond_apply(r);
        auto p = z.clone();

        T rs_old = z.dot(r);
        T y_norm = y.norm2() + (T)1e-10; // normalizer: ||y||_2
#ifdef DEBUG
        cudaProfilerStart();
#endif
        for (size_t iter = 0; iter < max_iter; iter++) {

            // Ap = A(p),  pAp = p^T Ap
            auto Ap = A(p);
            T pAp = Ap.dot(p);
            if (std::abs(pAp) < 1.e-10) {
                std::cerr << std::format(
                    "CG: p^T A p is too small ({:.5e}), stopping\n", pAp);
                break;
            }

            T alpha = rs_old / pAp;

            // step norm: ||delta_x|| = |alpha| * ||p|| (computed before updates)
            T dx = std::abs(alpha) * std::sqrt(p.dot(p)) /
                   (std::sqrt(x.dot(x)) + (T)1e-10);

            vec_xpay(x, p, alpha);   // x += alpha * p
            vec_xpay(r, Ap, -alpha); // r -= alpha * Ap
            // apply_support(x);

            // Apply preconditioner and compute new residual norm
            T rs_new = 0;
            z = precond_apply(r);
            rs_new = z.dot(r);

            if (std::abs(rs_old) < 2e-10) {
                std::cerr << "rs_old near zero, CG stagnated\n";
                break;
            }

            // Update search direction: p = z + beta * p
            T beta = rs_new / rs_old;
            vec_axpy(p, beta, z);
            rs_old = rs_new;

            T res = r.norm2() / y_norm;
            if (logger)
                logger->log(std::format(
                    "\tCG iter {:5d}: residual = {:.5e}, ||dx|| = {:.5e}\n",
                    iter + 1, res, dx));
            if (res < tol || dx < xtol) break;
        }
#ifdef DEBUG
        cudaProfilerStop();
#endif
        return x;
    }

    // Template instantiations
    template VecArray<float> cgsolver(const gpuFunction<float> &A,
                                      const VecArray<float> &y,
                                      const VecArray<float> &x0, size_t max_iter,
                                      float tol, float xtol, dims_t support_dims,
                                      float lambda, Logger *logger);
    template VecArray<double> cgsolver(const gpuFunction<double> &A,
                                       const VecArray<double> &y,
                                       const VecArray<double> &x0, size_t max_iter,
                                       double tol, double xtol, dims_t support_dims,
                                       double lambda, Logger *logger);

} // namespace tomocam::gpu::opt
