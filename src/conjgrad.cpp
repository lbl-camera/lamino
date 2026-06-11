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
#include "optimize.h"

namespace tomocam::opt {

    template <typename T>
    VecArray<T> cgsolver(const Function<T> &A, const VecArray<T> &y,
                         const VecArray<T> &x0, size_t max_iter, T tol, T xtol,
                         T lambda) {

        // initialize
        VecArray<T> x = x0.clone();
        VecArray<T> r = VecArray<T>::zeros(x0[0].dims());
        VecArray<T> p;

        // placeholder for preconditioner, currently identity
        auto precond_apply = [](const Array<T> &r) { return r.clone(); };

        VecArray<T> tmp = A(x);
        for (size_t i = 0; i < 3; i++) { r[i] = y[i] - tmp[i]; }

        VecArray<T> z;
        for (size_t i = 0; i < 3; i++) {
            z[i] = precond_apply(r[i]);
            p[i] = z[i].clone();
        }
        T rs_old = z.dot(r);

        for (size_t iter = 0; iter < max_iter; iter++) {

            // compute Ap
            VecArray<T> Ap = Ad(p);
            T pAp = p.dot(Ap);
            if (std::abs(pAp) < 1.e-10) {
                std::cerr << "pAp is close to zero\n";
                break;
            }
            T alpha = rs_old / pAp;

            // step norm: ||delta_x|| = |alpha| * ||p||
            T pp = 0;
            for (size_t i = 0; i < 3; i++) { pp += array::dot(p[i], p[i]); }
            T dx_norm = std::abs(alpha) * std::sqrt(pp);

            for (size_t i = 0; i < 3; i++) {
                x[i] += p[i] * alpha;
                r[i] -= Ap[i] * alpha;
            }

            // apply preconditioner
            for (size_t i = 0; i < 3; i++) { z[i] = precond.apply(r[i]); }
            T rs_new = z.dot(r);
            std::cout << std::format("\tCG iter: {}, residual: {}\n", iter,
                                     std::sqrt(rs_new));
            if (std::sqrt(rs_new) < tol) { break; }

            for (size_t i = 0; i < 3; i++) {
                p[i] = z[i] + p[i] * (rs_new / rs_old);
            }
            rs_old = rs_new;

#ifdef DEBUG
            // calculate ratio of data-fidelity to regularization
            T data_fidelity = 0;
            T regularization = 0;
            auto Atx = A(x);
            std::array<Array<T>, 3> x_std = {x[0].clone(), x[1].clone(),
                                             x[2].clone()};
            auto Htx = demag(x_std);
            for (size_t i = 0; i < 3; i++) {
                data_fidelity += array::norm2(Atx[i] - y[i]);
                regularization += array::norm2(Htx[i] * lambda);
            }
            std::cout << std::format(
                "\tCG iter: {:5}, residual: {:.5e}, dx: {:.5e}\n", iter,
                std::sqrt(res), dx_norm);
            if (std::sqrt(res) < tol || dx_norm < xtol) { break; }
        }
        return x;
    }

    // template instantiations
    template VecArray<float> cgsolver<float>(const Function<float> &A,
                                             const VecArray<float> &y,
                                             const VecArray<float> &x0,
                                             size_t max_iter, float tol, float xtol,
                                             float lambda);
    template VecArray<double> cgsolver<double>(const Function<double> &A,
                                               const VecArray<double> &y,
                                               const VecArray<double> &x0,
                                               size_t max_iter, double tol,
                                               double xtol, double lambda);

} // namespace tomocam::opt
