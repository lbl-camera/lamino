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

#include <algorithm>
#include <execution>
#include <format>

#include "array.h"
#include "array_ops.h"
#include "bregman.h"
#include "optimize.h"

namespace tomocam::opt {

    constexpr double EPSILON = 1e-8;

    /* Solve the LS problem with TV regularization:
     * min_x 0.5 ||A x - y||^2 + lambda TV(x)
     * using split Bregman method
     */
    template <typename T>
    VecArray<T> split_bregman(const std::vector<Function<T>> &As,
                              const std::vector<VecArray<T>> &yTs,
                              const VecArray<T> &x0, T lambda, T mu,
                              size_t outer_max, size_t inner_max, T tol, T xtol) {

        // Initialize variables
        auto dims = x0[0].dims();
        std::array<Array<T>, 3> x;
        std::array<Array<T>, 3> x_old;
        for (size_t i = 0; i < 3; ++i) {
            x[i] = x0[i].clone();
            x_old[i] = x0[i].clone();
        }

        // aux variable d and b are 3x3 matrices since gradient of vector is a 3x3
        // matrix
        std::array<std::array<Array<T>, 3>, 3> d;
        std::array<std::array<Array<T>, 3>, 3> b;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                d[i][j] = Array<T>::zeros(dims);
                b[i][j] = Array<T>::zeros(dims);
            }
        }

        // update A^TA to add laplacian of x
        // Ap  = (A^TA  +   ∇^T∇) u
        Function<T> Ap = [&](const std::array<Array<T>, 3> &u) {
            std::array<Array<T>, 3> result;
            for (size_t i = 0; i < 3; ++i) { result[i] = laplacian(u[i]) * mu; }
            for (auto &A : As) { result += A(u); }
            return result;
        };

        for (size_t iter = 0; iter < outer_max; ++iter) {

            // update RHS := R^T y1 + R^Ty2 + μ∇^T(d - b)
            std::array<Array<T>, 3> rhs;
            for (size_t i = 0; i < 3; ++i) { rhs[i] = divergence(d[i] - b[i]) * mu; }
            for (auto &yT : yTs) {
                for (size_t i = 0; i < 3; ++i) { rhs[i] += yT[i]; }
            }

            // use conjugate gradient to solve the linear system
            x = cgsolver<T>(Ap, rhs, x_old, inner_max, tol);

            /* Isotropic TV shrinkage */
            // compute gradient of solution
            std::array<std::array<Array<T>, 3>, 3> grad_x;
            for (size_t i = 0; i < 3; ++i) { grad_x[i] = grad_u(x[i]); }

            std::array<Array<T>, 3> sk;
            // compute shrinkage factor
            for (size_t i = 0; i < 3; ++i) {
                sk[i] = Array<T>::zeros(dims);
                sk[i] += EPSILON; // to avoid division by zero
                for (size_t j = 0; j < 3; ++j) {
                    sk[i] += (grad_x[i][j] + b[i][j]) * (grad_x[i][j] + b[i][j]);
                }
                std::transform(std::execution::par_unseq, sk[i].begin(), sk[i].end(),
                               sk[i].begin(), [](T v) { return std::sqrt(v); });
            }
            // update d with shrinkage
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    auto temp = (grad_x[i][j] + b[i][j]) / sk[i];
                    std::transform(std::execution::par_unseq, temp.begin(),
                                   temp.end(), sk[i].begin(), d[i][j].begin(),
                                   [lambda, mu](T val, T sk_val) {
                                       return std::max(sk_val - lambda / mu, T(0)) *
                                              val;
                                   });
                }
            }

            // Bregman update
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) { b[i][j] += grad_x[i][j] - d[i][j]; }
            }

            // Check convergence
            T norm_diff = T(0);
            for (size_t i = 0; i < 3; ++i) {
                norm_diff += array::norm2(x[i] - x_old[i]) /
                             (array::norm2(x_old[i]) + static_cast<T>(EPSILON));
            }
            std::cout << std::format(
                "Outer iter: {}, ‖xᵏ⁺¹ − xᵏ‖₂ / ‖xᵏ‖₂: {:.6e}\n", iter, norm_diff);
            for (size_t i = 0; i < 3; ++i) { x_old[i] = std::move(x[i]); }
            if (norm_diff < xtol) { break; }
        }
        return x;
    }
    // Explicit template instantiation for float and double
    template VecArray<float> split_bregman(const std::vector<Function<float>> &As,
                                           const std::vector<VecArray<float>> &yTs,
                                           const VecArray<float> &x0, float lambda,
                                           float mu, size_t outer_max,
                                           size_t inner_max, float tol, float xtol);
    template VecArray<double> split_bregman(const std::vector<Function<double>> &As,
                                            const std::vector<VecArray<double>> &yTs,
                                            const VecArray<double> &x0,
                                            double lambda, double mu,
                                            size_t outer_max, size_t inner_max,
                                            double tol, double xtol);

} // namespace tomocam::opt
