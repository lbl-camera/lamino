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

#include "array.h"
#include "optimize.h"

namespace tomocam::opt {

    /**! Conjugate Gradient Least Squares (CGLS) algorithm for solving the
     * linear system /f$ [A_1, A_2]^T x = [b_1, b_2]^T /f$.
     * \param A1 Function representing the first linear operator
     * \param A2 Function representing the second linear operator
     * \param b1 Right-hand side vector corresponding to A1
     * \param b2 Right-hand side vector corresponding to A2
     * \param x0 Initial guess for the solution
     * \param max_iters Maximum number of iterations
     * \param tol Tolerance for convergence
     * \return Approximate solution vector x
     */
    template <typename T>
    VecArray<T> cgls(const Function<T> &A1, const Function<T> &A2,
                     const VecArray<T> &b1, const VecArray<T> &b2,
                     const VecArray<T> &x0, const size_t max_iters, const T tol) {

        // initialize
        VecArray<T> x;
        for (size_t i = 0; i < x0.size(); i++) { x[i] = x0[i].clone(); }

        VecArray<T> r = (b1 - A1(x)) + (b2 - A2(x));
        VecArray<T> p;
        for (size_t i = 0; i < r.size(); i++) { p[i] = r[i].clone(); }

        T r_norm = (T)0;
        for (size_t i = 0; i < r.size(); i++) { r_norm += array::dot(r[i], r[i]); }

        // iterate
        for (size_t iter = 0; iter < max_iters; iter++) {

            // compute Ap
            VecArray<T> Ap = A1(p) + A2(p);

            // compute alpha
            T pAp = (T)0;
            for (size_t i = 0; i < r.size(); i++) { pAp += array::dot(p[i], Ap[i]); }
            T alpha = r_norm / pAp;

            // update x
            VecArray<T> r_new;
            for (size_t i = 0; i < x.size(); i++) {
                x[i] += p[i] * alpha;
                r_new[i] = r[i] - Ap[i] * alpha;
            }

            // check convergence
            T res_norm = (T)0;
            for (size_t i = 0; i < r_new.size(); i++) {
                res_norm += array::dot(r_new[i], r_new[i]);
            }
            if (res_norm < tol * tol) { break; }

            // update r and p
            T beta = res_norm / r_norm;
            for (size_t i = 0; i < p.size(); i++) {
                p[i] = r_new[i] + p[i] * beta;
                r[i] = std::move(r_new[i]);
            }
            r_norm = res_norm;
        }

        return x;
    }

    // explicit template instantiations
    template VecArray<float>
    cgls<float>(const Function<float> &A1, const Function<float> &A2,
                const VecArray<float> &b1, const VecArray<float> &b2,
                const VecArray<float> &x0, const size_t max_iters, const float tol);
    template VecArray<double> cgls<double>(const Function<double> &A1,
                                           const Function<double> &A2,
                                           const VecArray<double> &b1,
                                           const VecArray<double> &b2,
                                           const VecArray<double> &x0,
                                           const size_t max_iters, const double tol);
} // namespace tomocam::opt
