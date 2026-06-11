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

#include <concepts>
#include <format>

#include <thrust/iterator/zip_iterator.h>
#include <thrust/transform.h>

#include "gpu/bregman.h"
#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/vec_array.h"
#include "gpu/finitediff.h"
#include "gpu/gpu_opt.h"
#include "gpu/mem_check.h"

namespace tomocam::gpu::opt {

    constexpr double EPSILON = 1e-8;

    // isotropic TV shrinkage functor
    template <typename T>
    struct ShrinkFn {
        T lambda;
        T mu;
        ShrinkFn(T lambda, T mu) : lambda(lambda), mu(mu) {}

        template <typename Tuple>
        __host__ __device__ T operator()(const Tuple &t) const {
            using namespace cuda::std;
            T dx = thrust::get<0>(t);
            T b = thrust::get<1>(t);
            T sk = thrust::get<2>(t);
            T val = (dx + b) / (sk + (T)EPSILON);
            return max(T(0), sk - lambda / mu) * val;
        }
    };

    template <typename T>
    DeviceArray<T> shrink(const DeviceArray<T> &dx, const DeviceArray<T> &b,
                          const DeviceArray<T> &sk, T lambda, T mu) {
        DeviceArray<T> d(b.dims());
        auto zip = thrust::zip_iterator(
            thrust::make_tuple(dx.begin(), b.begin(), sk.begin()));
        thrust::transform(zip, zip + dx.size(), d.begin(), ShrinkFn(lambda, mu));
        return d;
    }

    /// compute sk = sqrt(∑(dx[i] + b[i])^2) for i=1,2,3
    template <typename T>
    struct AddSquareFn {
        __host__ __device__ T operator()(T dx, T b) const {
            return (dx + b) * (dx + b);
        }
    };
    template <typename T>
    DeviceArray<T> compute_sk(const std::array<VecArray<T>, 3> &dx,
                              const std::array<VecArray<T>, 3> &b) {
        auto dims = dx[0][0].dims();
        DeviceArray<T> sk(dims);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                DeviceArray<T> temp(dims);
                thrust::transform(dx[i][j].begin(), dx[i][j].end(), b[i][j].begin(),
                                  temp.begin(), AddSquareFn<T>());
                sk += temp;
            }
        }
        return array::sqrt(sk);
    }

    /// Split Bregman method for TV-regularized least squares problems
    template <typename T>
    VecArray<T> split_bregman(const gpuFunction<T> &A, const VecArray<T> &yT,
                              const VecArray<T> &x0, T lambda, T mu,
                              size_t outer_max, size_t inner_max, T tol, T xtol) {

        // initialize variables
        VecArray<T> x = x0.clone();
        VecArray<T> x_old = x0.clone();

        // bregman variables
        std::array<VecArray<T>, 3> b;
        std::array<VecArray<T>, 3> d;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                b[i][j] = DeviceArray<T>(x0[0].dims());
                d[i][j] = DeviceArray<T>(x0[0].dims());
            }
        }

        // update A^TA to add laplacian of x
        // Ap  = (A^TA + μ∇^T∇) u
        gpuFunction<T> Ap = [&A, &mu](const VecArray<T> &u) {
            auto Au = A(u);
            for (size_t i = 0; i < 3; ++i) {
                auto lap_u = neg_laplacian(u[i]);
                array::xpay(Au[i], lap_u, mu);
            }
            return Au;
        };

        for (int iter = 0; iter < outer_max; ++iter) {

            // x-update: solve (A^TA + μ∇^T∇)x = A^T y + μ∇^T(d - b)
            std::array<VecArray<T>, 3> d_b;
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) { d_b[i][j] = d[i][j] - b[i][j]; }
            }

            // compute right-hand side: A^T y + μ∇^T(d - b)
            auto rhs = yT.clone();
            for (size_t i = 0; i < 3; ++i) {
                auto neg_div_db = neg_divergence(d_b[i]);
                array::xpay(rhs[i], neg_div_db, mu);
            }

            // use conjugate gradient to solve the linear system
            x = cgsolver(Ap, rhs, x, inner_max, tol, xtol);

            // isotropic TV shrinkage
            std::array<VecArray<T>, 3> dx;
            for (size_t i = 0; i < 3; ++i) { dx[i] = std::move(grad(x[i])); }

            // sk = sqrt(∑(dx[i] + b[i])^2)
            auto sk = compute_sk(dx, b);

            // update d inplace
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    d[i][j] = shrink(dx[i][j], b[i][j], sk, lambda, mu);
                }
            }

            // Bregman update
            // b[i] = b[i] + (dx[i] - d[i])
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) { b[i][j] += dx[i][j] - d[i][j]; }
            }

            // Check convergence
            // norm_diff = ‖xᵏ⁺¹ − xᵏ‖₂ / ‖xᵏ‖₂
            T norm_diff = (x - x_old).norm2() / (x_old.norm2() + (T)EPSILON);
            std::cout << std::format(
                "Outer iter: {}, ‖xᵏ⁺¹ − xᵏ‖₂ / ‖xᵏ‖₂: {:.6e}\n", iter, norm_diff);
            x_old = x.clone();
            if (norm_diff < xtol) { break; }
        }
        return x;
    }

    // explicit template instantiations for float and double
    template VecArray<float> split_bregman(const gpuFunction<float> &A,
                                           const VecArray<float> &yT,
                                           const VecArray<float> &x0, float lambda,
                                           float mu, size_t outer_max,
                                           size_t inner_max, float tol, float xtol);
    template VecArray<double>
    split_bregman(const gpuFunction<double> &A, const VecArray<double> &yT,
                  const VecArray<double> &x0, double lambda, double mu,
                  size_t outer_max, size_t inner_max, double tol, double xtol);
#ifdef DEBUG
    // compute_sk test
    template DeviceArray<float>
    compute_sk(const std::array<VecArray<float>, 3> &dx,
               const std::array<VecArray<float>, 3> &b);
    template DeviceArray<double>
    compute_sk(const std::array<VecArray<double>, 3> &dx,
               const std::array<VecArray<double>, 3> &b);

    // shrink test
    template DeviceArray<float> shrink(const DeviceArray<float> &dx,
                                       const DeviceArray<float> &b,
                                       const DeviceArray<float> &sk, float lambda,
                                       float mu);
    template DeviceArray<double> shrink(const DeviceArray<double> &dx,
                                        const DeviceArray<double> &b,
                                        const DeviceArray<double> &sk, double lambda,
                                        double mu);
#endif
} // namespace tomocam::gpu::opt
