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

#include <array>
#include <execution>
#include <format>
#include <functional>

#include "array.h"
#include "array_ops.h"
#include "bregman.h"
#include "mask.h"
#include "optimize.h"

namespace tomocam::opt {

    template <typename T>
    VecArray<T> cgsolver(const Function<T> &A, const VecArray<T> &y,
                         const VecArray<T> &x0, size_t max_iter, T tol, T xtol,
                         const Array<T> &sup_mask, T lambda, Logger *logger) {

        // initialize
        VecArray<T> x = clone(x0);

        // placeholder for preconditioner, currently identity
        auto precond_apply = [](const VecArray<T> &r) { return clone(r); };

        // support constraint
        auto apply_support = [&sup_mask](VecArray<T> &v) {
            for (size_t i = 0; i < 3; i++) { v[i] *= sup_mask; }
        };

        // project x0 onto support subspace once so every iterate stays in it
        apply_support(x);

        T y_norm = std::sqrt(dot(y, y)) + (T)1e-10;

        // compute initial residual
        VecArray<T> r = y - A(x);
        VecArray<T> z = precond_apply(r);
        VecArray<T> p = clone(z);
        T rs_old = dot(z, r);

        for (size_t iter = 0; iter < max_iter; iter++) {

            VecArray<T> Ap = A(p);
            T pAp = dot(p, Ap);
            if (std::abs(pAp) < 2.e-10) {
                std::cerr << "pAp is close to zero\n";
                break;
            }
            T alpha = rs_old / pAp;

            // step norm: ||delta_x|| = |alpha| * ||p||
            T pp = dot(p, p);
            T xnorm = dot(x, x);
            T dx = std::abs(alpha) * std::sqrt(pp) / (std::sqrt(xnorm) + (T)1.e-10);

            for (size_t i = 0; i < 3; i++) {
                x[i] += p[i] * alpha;
                r[i] -= Ap[i] * alpha;
            }
            apply_support(x);

            // apply preconditioner
            z = precond_apply(r);
            T rs_new = dot(z, r);

            // update p
            if (std::abs(rs_old) < 2.e-10) {
                std::cerr << "rs_old near zero, CG stagnated\n";
                break;
            }
            T beta = rs_new / rs_old;
            for (size_t i = 0; i < 3; i++) { p[i] = z[i] + p[i] * beta; }
            rs_old = rs_new;

            T res = std::sqrt(dot(r, r)) / y_norm;
            if (logger)
                logger->log(std::format(
                    "\tCG iter: {:5}, residual: {:.5e}, ||dx||: {:.5e}\n", iter + 1,
                    res, dx));
            if (res < tol || dx < xtol) { break; }
        }
        return x;
    }

    // template instantiations
    template VecArray<float> cgsolver<float>(const Function<float> &A,
                                             const VecArray<float> &y,
                                             const VecArray<float> &x0,
                                             size_t max_iter, float tol, float xtol,
                                             const Array<float> &sup_mask,
                                             float lambda, Logger *logger);
    template VecArray<double>
    cgsolver<double>(const Function<double> &A, const VecArray<double> &y,
                     const VecArray<double> &x0, size_t max_iter, double tol,
                     double xtol, const Array<double> &sup_mask, double lambda,
                     Logger *logger);

} // namespace tomocam::opt
