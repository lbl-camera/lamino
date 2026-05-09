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
#include <iostream>

#include "array.h"
#include "array_ops.h"
#include "bregman.h"
#include "demag.h"
#include "optimize.h"
#include "support.h"

namespace tomocam::opt {

    template <typename T>
    VecArray<T> cgsolver(const Function<T> &A, const VecArray<T> &y,
                         const VecArray<T> &x0, size_t max_iter, T tol, T lambda) {

        auto make_vecarray = [](const dims_t &dims) {
            return VecArray<T>{Array<T>(dims), Array<T>(dims), Array<T>(dims)};
        };
        // initialize
        VecArray<T> x = {x0[0].clone(), x0[1].clone(), x0[2].clone()};
        auto r = make_vecarray(x0[0].dims());
        auto p = make_vecarray(x0[0].dims());
        auto z = make_vecarray(x0[0].dims());

        // TODO: apply lambda regularization (A + lambda*C)x
        Function<T> Ad = [&A, lambda](const VecArray<T> &x) {
            VecArray<T> Ax = A(x);
            return Ax;
        };

        // placeholder preconditioner (identity)
        auto precond_apply = [](const Array<T> &r) { return r.clone(); };

        VecArray<T> tmp = Ad(x);
        for (size_t i = 0; i < 3; i++) { r[i] = y[i] - tmp[i]; }

        T rs_old = 0;
        for (size_t i = 0; i < 3; i++) {
            z[i] = precond_apply(r[i]);
            p[i] = z[i].clone();
            rs_old += array::dot(z[i], r[i]);
        }

        for (size_t iter = 0; iter < max_iter; iter++) {

            // project search direction before applying operator
            VecArray<T> Ap = Ad(p);
            T pAp = 0;
            for (size_t i = 0; i < 3; i++) { pAp += array::dot(p[i], Ap[i]); }
            if (std::abs(pAp) < 1.e-10) {
                std::cerr << "pAp is close to zero\n";
                break;
            }

            T alpha = rs_old / pAp;
            for (size_t i = 0; i < 3; i++) {
                x[i] += p[i] * alpha;
                r[i] -= Ap[i] * alpha;
            }

            // apply preconditioner
            T rs_new = 0;
            for (size_t i = 0; i < 3; i++) {
                z[i] = precond_apply(r[i]);
                rs_new += array::dot(z[i], r[i]);
            }

            for (size_t i = 0; i < 3; i++) {
                p[i] = z[i] + p[i] * (rs_new / rs_old);
            }
            rs_old = rs_new;

            T res = 0;
            for (size_t i = 0; i < 3; i++) { res += array::dot(r[i], r[i]); }
            std::cout << std::format("\tCG iter: {}, residual: {}\n", iter,
                                     std::sqrt(res));
            if (std::sqrt(res) < tol) { break; }
        }
        return x;
    }

    // template instantiations
    template VecArray<float> cgsolver<float>(const Function<float> &A,
                                             const VecArray<float> &y,
                                             const VecArray<float> &x0,
                                             size_t max_iter, float tol,
                                             float lambda);
    template VecArray<double> cgsolver<double>(const Function<double> &A,
                                               const VecArray<double> &y,
                                               const VecArray<double> &x0,
                                               size_t max_iter, double tol,
                                               double lambda);

} // namespace tomocam::opt
