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

#include "array.h"
#include "array_ops.h"
#include "tomocam.h"

/* TODO:
 * implement finite difference for Bregman TV
 * implment CG solver for the linear system
 */
namespace tomocam::opt {

    template <typename T>
    Array<T> shrink(const Array<T> &v, T kappa) {
        Array<T> result = v.clone();
        for (size_t i = 0; i < v.size(); ++i) {
            T val = v[i];
            if (val > kappa) {
                result[i] = val - kappa;
            } else if (val < -kappa) {
                result[i] = val + kappa;
            } else {
                result[i] = 0;
            }
        }
        return result;
    }

    template <typename T>
    Array<T> split_bregman(std::function<Array<T>(const Array<T> &)> A,
                           const Array<T> &yT, T lambda, T mu, size_t outer_max,
                           size_t inner_max, T tol) {

        Array<T> x = yT.clone();
        Array<T> d = Array<T>::zeros(x.dims());
        Array<T> b = Array<T>::zeros(x.dims());

        for (int iter = 0; iter < outer_max; ++iter) {

            // x-update: solve (A^T A + mu I)x = A^T b + mu (d - bregman)
            auto rhs = yT + mu * (d - b);

            // user conjugate gradient to solve the linear system
            auto x_cg = cgsolve(A, rhs, mu, inner_max, tol);

            // d-update: shrinkage step
            d = shrink(x_cg + b, lambda / mu);

            // Bregman update
            b += (x - d);

            // Check convergence
            T norm_diff = array::norm2(x - d);
            if (norm_diff < tol) { break; }
        }
        return x;
    }
} // namespace tomocam::opt
