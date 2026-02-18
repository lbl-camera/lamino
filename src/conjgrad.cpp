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
#include "optimize.h"

namespace tomocam::opt {

    template <typename T>
    VecArray<T> cgsolver(const Function<T> &A, const VecArray<T> &y,
                         const VecArray<T> &x0, size_t max_iter, T tol) {

        // initialize
        VecArray<T> x;
        VecArray<T> r;
        VecArray<T> p;
        for (size_t i = 0; i < 3; i++) { x[i] = x0[i].clone(); }
        VecArray<T> tmp = A(x);
        for (size_t i = 0; i < 3; i++) {
            r[i] = y[i] - tmp[i];
            p[i] = r[i].clone();
        }
        T rs_old = 0;
        for (size_t i = 0; i < 3; i++) { rs_old += array::dot(r[i], r[i]); }

        for (size_t iter = 0; iter < max_iter; iter++) {

            VecArray<T> Ap = A(p);
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

            T rs_new = 0;
            for (size_t i = 0; i < 3; i++) { rs_new += array::dot(r[i], r[i]); }
            std::cout << std::format("iter: {}, residual: {}\n", iter,
                                     std::sqrt(rs_new));
            if (std::sqrt(rs_new) < tol) { break; }

            for (size_t i = 0; i < 3; i++) {
                p[i] = r[i] + p[i] * (rs_new / rs_old);
            }
            rs_old = rs_new;
        }
        return x;
    }

    // template instantiations
    template VecArray<float> cgsolver<float>(const Function<float> &A,
                                             const VecArray<float> &y,
                                             const VecArray<float> &x0,
                                             size_t max_iter, float tol);
    template VecArray<double> cgsolver<double>(const Function<double> &A,
                                               const VecArray<double> &y,
                                               const VecArray<double> &x0,
                                               size_t max_iter, double tol);

} // namespace tomocam::opt
