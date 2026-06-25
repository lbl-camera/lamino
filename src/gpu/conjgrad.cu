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
#include <iostream>
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

    // -------------------------------------------------------------------------
    // GPU Conjugate Gradient Solver
    // Same preconditioned CG algorithm as src/conjgrad.cpp, with GPU-native
    // DeviceArray types and Thrust/cuFFT inner operations.
    // -------------------------------------------------------------------------
    template <typename T>
    VecArray<T> cgsolver(const gpuFunction<T> &A, const VecArray<T> &y,
                         const VecArray<T> &x0, size_t max_iter, T tol, T xtol,
                         dims_t support_dims, T lambda) {

        auto precond_apply = [](const DeviceArray<T> &r) { return r.clone(); };

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
        VecArray<T> z{precond_apply(r[0]), precond_apply(r[1]), precond_apply(r[2])};
        auto p = z.clone();

        T rs_old = z.dot(r);
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
            vec_xpay(x, p, alpha);   // x += alpha * p
            apply_support(x);
            vec_xpay(r, Ap, -alpha); // r -= alpha * Ap

            // Apply preconditioner and compute new residual norm
            T rs_new = 0;
            for (size_t i = 0; i < 3; ++i) { z[i] = precond_apply(r[i]); }
            rs_new = z.dot(r);

            // Update search direction: p = z + beta * p
            T beta = rs_new / rs_old;
            vec_axpy(p, beta, z);
            rs_old = rs_new;

            T res = r.norm2();
            if (res < tol) break;

            // dx: relative step size (computed before p is updated)
            T dx = std::abs(alpha) * std::sqrt(p.dot(p)) /
                   (std::sqrt(x.dot(x)) + (T)1e-10);

            if (dx < xtol) {
                std::cout << "CG converged based on solution change\n";
                break;
            }
            std::cout << std::format(
                "\t CG iter {:3d}: residual = {:.6e}, dx = {:.6e}\n", iter + 1, res,
                dx);
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
                                      float lambda);
    template VecArray<double> cgsolver(const gpuFunction<double> &A,
                                       const VecArray<double> &y,
                                       const VecArray<double> &x0, size_t max_iter,
                                       double tol, double xtol, dims_t support_dims,
                                       double lambda);

} // namespace tomocam::gpu::opt
