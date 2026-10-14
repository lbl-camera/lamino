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

#include <cmath>
#include <format>
#include <stdexcept>

#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/transform_reduce.h>

#include "gpu/bregman.h"
#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/finitediff.h"
#include "gpu/gpu_opt.h"
#include "gpu/vec_array.h"

namespace tomocam::gpu::opt {

    constexpr double MAG_EPSILON = 1e-12;

    namespace {
        // Binary support mask on the device (mask > 0.5 inside). Voxels outside
        // it are padding and are left unconstrained.
        template <typename T>
        struct Support {
            const T *m;   // device pointer to the mask
            size_t count; // number of voxels inside

            __host__ __device__ bool contains(size_t i) const {
                return m[i] > static_cast<T>(0.5);
            }
        };

        // copy of x0, zeroed outside the support if zero_outside
        template <typename T>
        VecArray<T> initial_guess(const VecArray<T> &x0, const DeviceArray<T> &mask,
                                  bool zero_outside) {
            VecArray<T> x = x0.clone();
            if (zero_outside) {
                for (size_t i = 0; i < 3; ++i) { x[i] *= mask; }
            }
            return x;
        }

        template <typename T>
        struct InsideFn {
            __host__ __device__ bool operator()(T v) const {
                return v > static_cast<T>(0.5);
            }
        };

        template <typename T>
        Support<T> make_support(const DeviceArray<T> &mask, const dims_t &dims) {
            if (!(mask.dims() == dims)) {
                throw std::invalid_argument(
                    "constant-magnitude mask must have the volume dimensions");
            }
            const T *m = thrust::raw_pointer_cast(mask.begin());
            size_t count = thrust::count_if(thrust::device, mask.begin(), mask.end(),
                                            InsideFn<T>{});
            if (count == 0) {
                throw std::invalid_argument("constant-magnitude mask is empty");
            }
            return Support<T>{m, count};
        }

        // |v| inside the mask, 0 outside
        template <typename T>
        struct MagnitudeFn {
            Support<T> roi;

            template <typename Tuple>
            __host__ __device__ double operator()(const Tuple &t) const {
                if (!roi.contains(thrust::get<0>(t))) return 0.0;
                double x = thrust::get<1>(t);
                double y = thrust::get<2>(t);
                double z = thrust::get<3>(t);
                return sqrt(x * x + y * y + z * z);
            }
        };

        // rescale v to |v| = M0 inside the mask. A voxel with |v| <= eps keeps
        // the previous value q (if given), otherwise becomes (dflt_x, 0, 0).
        // Outside the mask v is kept, or set to 0 if zero_outside.
        template <typename T>
        struct ProjectFn {
            Support<T> roi;
            double M0;
            const T *qx, *qy, *qz;
            double dflt_x;
            bool zero_outside;
            double eps;

            template <typename Tuple>
            __host__ __device__ void operator()(Tuple t) const {
                size_t i = thrust::get<0>(t);
                if (!roi.contains(i)) {
                    if (zero_outside) {
                        thrust::get<1>(t) = T(0);
                        thrust::get<2>(t) = T(0);
                        thrust::get<3>(t) = T(0);
                    }
                    return;
                }
                double x = thrust::get<1>(t);
                double y = thrust::get<2>(t);
                double z = thrust::get<3>(t);
                double mag = sqrt(x * x + y * y + z * z);
                if (mag > eps) {
                    double s = M0 / mag;
                    thrust::get<1>(t) = static_cast<T>(x * s);
                    thrust::get<2>(t) = static_cast<T>(y * s);
                    thrust::get<3>(t) = static_cast<T>(z * s);
                } else if (qx) {
                    thrust::get<1>(t) = qx[i];
                    thrust::get<2>(t) = qy[i];
                    thrust::get<3>(t) = qz[i];
                } else {
                    thrust::get<1>(t) = static_cast<T>(dflt_x);
                    thrust::get<2>(t) = T(0);
                    thrust::get<3>(t) = T(0);
                }
            }
        };

        template <typename T>
        T mean_magnitude(const VecArray<T> &v, const Support<T> &roi) {
            auto zip = thrust::make_zip_iterator(
                thrust::make_tuple(thrust::counting_iterator<size_t>(0),
                                   v[0].begin(), v[1].begin(), v[2].begin()));
            double sum =
                thrust::transform_reduce(zip, zip + v[0].size(), MagnitudeFn<T>{roi},
                                         0.0, thrust::plus<double>());
            return static_cast<T>(sum / static_cast<double>(roi.count));
        }

        template <typename T>
        void project_magnitude_inplace(VecArray<T> &v, T M0, const Support<T> &roi,
                                       const VecArray<T> *z_prev = nullptr,
                                       double dflt_x = 0.0,
                                       bool zero_outside = false,
                                       double eps = MAG_EPSILON) {
            auto zip = thrust::make_zip_iterator(
                thrust::make_tuple(thrust::counting_iterator<size_t>(0),
                                   v[0].begin(), v[1].begin(), v[2].begin()));
            auto raw = [&](size_t c) -> const T * {
                return z_prev ? thrust::raw_pointer_cast((*z_prev)[c].begin())
                              : nullptr;
            };
            thrust::for_each(zip, zip + v[0].size(),
                             ProjectFn<T>{roi, static_cast<double>(M0), raw(0),
                                          raw(1), raw(2), dflt_x, zero_outside,
                                          eps});
        }
    } // namespace

    template <typename T>
    T mean_magnitude(const VecArray<T> &v) {
        auto ones = DeviceArray<T>::ones(v[0].dims());
        return mean_magnitude(v, make_support(ones, v[0].dims()));
    }

    template <typename T>
    VecArray<T> project_magnitude(const VecArray<T> &v, T M0) {
        VecArray<T> rv = v.clone();
        auto ones = DeviceArray<T>::ones(v[0].dims());
        project_magnitude_inplace(rv, M0, make_support(ones, v[0].dims()));
        return rv;
    }

    /// Constant magnitude + L2 smoothing, mirrors src/regularized_admm.cpp
    template <typename T>
    VecArray<T> project_sphere(const VecArray<T> &v, T Ms, const VecArray<T> &z_prev,
                               T eps) {
        VecArray<T> rv = v.clone();
        auto ones = DeviceArray<T>::ones(v[0].dims());
        project_magnitude_inplace(rv, Ms, make_support(ones, v[0].dims()), &z_prev,
                                  0.0, false, static_cast<double>(eps));
        return rv;
    }

    template <typename T>
    VecArray<T> constmag_l2(const gpuFunction<T> &A, const VecArray<T> &yT,
                            const VecArray<T> &x0, T lambda_l2, T mu,
                            size_t outer_max, size_t inner_max, T tol, T xtol,
                            const DeviceArray<T> &mask, Logger *logger,
                            bool zero_outside) {

        auto dims = x0[0].dims();
        if (!(dims == yT[0].dims())) {
            throw std::invalid_argument(
                "yT is not the raw data, but backprojected data, so it should have "
                "the same dimensions as x0");
        }
        const auto roi = make_support(mask, dims);
        const DeviceArray<T> no_support;

        // the magnitude is uniform but unknown: start from the mean magnitude
        // of the initial guess and re-estimate it in every z-update
        T M0 = mean_magnitude(x0, roi);
        if (logger)
            logger->log(std::format("Initial constant magnitude M0: {:.6e}\n", M0));

        VecArray<T> x;
        const VecArray<T> x_init = initial_guess(x0, mask, zero_outside);
        VecArray<T> x_old = x_init.clone();
        VecArray<T> z = x_init.clone();
        project_magnitude_inplace(z, M0, roi,
                                  static_cast<const VecArray<T> *>(nullptr), 0.0,
                                  zero_outside);
        VecArray<T> b{DeviceArray<T>(dims), DeviceArray<T>(dims),
                      DeviceArray<T>(dims)};

        // Ap = (A^TA + lambda_l2 ∇^T∇ + mu I) u
        gpuFunction<T> Ap = [&A, lambda_l2, mu](const VecArray<T> &u) {
            auto Au = A(u);
            for (size_t i = 0; i < 3; ++i) {
                auto lap_u = neg_laplacian(u[i]);
                array::xpay(Au[i], lap_u, lambda_l2);
                array::xpay(Au[i], u[i], mu);
            }
            return Au;
        };

        for (size_t iter = 0; iter < outer_max; ++iter) {

            {
                // rhs := R^T y + mu (z - b)
                auto rhs = yT.clone();
                vec_xpay(rhs, z, mu);
                vec_xpay(rhs, b, -mu);

                // x-update
                x = cgsolver(Ap, rhs, x_old, inner_max, tol, tol,
                             zero_outside ? mask : no_support, (T)0, nullptr);
            }

            // z-update: project (x + b) onto {|z| = M0 inside the mask, M0 free}.
            // The optimal M0 is the mean magnitude of (x + b) over the mask;
            // outside the mask z = x + b, so b stays zero there and the padding
            // is left unconstrained.
            z = x.clone();
            z += b;
            M0 = mean_magnitude(z, roi);
            project_magnitude_inplace(z, M0, roi,
                                      static_cast<const VecArray<T> *>(nullptr), 0.0,
                                      zero_outside);

            // dual update
            b += x;
            b -= z;

            // convergence: relative change of x and primal residual ||x - z||
            T xnorm = x.norm2() + (T)MAG_EPSILON;
            T rel_dx = (x - x_old).norm2() / (x_old.norm2() + (T)MAG_EPSILON);
            T rel_r = (x - z).norm2() / xnorm;
            if (logger)
                logger->log(std::format("Outer iter: {}, M0: {:.6e}, "
                                        "||dx||/||x||: {:.6e}, "
                                        "||x-z||/||x||: {:.6e}\n",
                                        iter, M0, rel_dx, rel_r));
            x_old = std::move(x);
            if (rel_dx < xtol && rel_r < xtol) { break; }
        }
        // return the feasible variable (as on the CPU)
        return z;
    }

    /// Combined constant magnitude / TV / L2, mirrors src/regularized_admm.cpp
    template <typename T>
    VecArray<T> regularized_admm(const gpuFunction<T> &A, const VecArray<T> &yT,
                                 const VecArray<T> &x0, const RegWeights<T> &w, T mu,
                                 size_t outer_max, size_t inner_max, T tol, T xtol,
                                 const DeviceArray<T> &mask, Logger *logger,
                                 const MagOptions<T> &mag, AdmmHistory<T> *history) {

        auto dims = x0[0].dims();
        if (!(dims == yT[0].dims())) {
            throw std::invalid_argument(
                "yT is not the raw data, but backprojected data, so it should have "
                "the same dimensions as x0");
        }
        if (!w.const_mag && !w.tv() && !w.l2()) {
            throw std::invalid_argument(
                "regularized_admm: enable at least one of const_mag, TV or L2");
        }
        if (mag.rho < T(0) || mag.rho_growth < T(1) || mag.rho_max < T(0)) {
            throw std::invalid_argument("regularized_admm: rho must be >= 0, "
                                        "rho_growth >= 1 (rho only increases)");
        }
        const bool use_mag = w.const_mag;
        const bool use_tv = w.tv();
        const auto roi = make_support(mask, dims);
        const DeviceArray<T> no_support;

        VecArray<T> x;
        const VecArray<T> x_init = initial_guess(x0, mask, mag.zero_outside);
        VecArray<T> x_old = x_init.clone();

        // magnitude splitting: z = x, dual b, penalty rho (continued upwards)
        T rho = mag.rho > T(0) ? mag.rho : mu;
        T M0 = 0;
        VecArray<T> z, b;
        if (use_mag) {
            M0 = (!mag.estimate_Ms && mag.Ms > T(0)) ? mag.Ms
                                                     : mean_magnitude(x0, roi);
            if (logger)
                logger->log(std::format("Initial constant magnitude Ms: {:.6e}, "
                                        "rho: {:.6e}\n",
                                        M0, rho));
            // z0 = P(x0); zero vectors in the mask get the default direction +x
            z = x_init.clone();
            project_magnitude_inplace(z, M0, roi,
                                      static_cast<const VecArray<T> *>(nullptr),
                                      static_cast<double>(M0), mag.zero_outside);
            b = VecArray<T>{DeviceArray<T>(dims), DeviceArray<T>(dims),
                            DeviceArray<T>(dims)};
            if (mag.zero_outside) {
                // z0 is 0 outside the mask, so start with the matching multiplier
                b = x_init.clone();
                b -= z;
            }
        }

        // TV splitting: d = grad x, dual bt (3x3 volumes)
        std::array<VecArray<T>, 3> d, bt;
        if (use_tv) {
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    d[i][j] = DeviceArray<T>(dims);
                    bt[i][j] = DeviceArray<T>(dims);
                }
            }
        }

        // Ap = (A^TA + (lambda_l2 + mu [tv]) ∇^T∇ + rho [mag] I) u
        const T lap_weight = w.lambda_l2 + (use_tv ? mu : T(0));
        const bool use_lap = use_tv || w.l2();
        // beta only enters the system matrix, never the right-hand side. rho is
        // captured by reference: it may grow between outer iterations.
        gpuFunction<T> Ap = [&A, lap_weight, use_lap, use_mag,
                             &rho](const VecArray<T> &u) {
            auto Au = A(u);
            for (size_t i = 0; i < 3; ++i) {
                if (use_lap) {
                    auto lap_u = neg_laplacian(u[i]);
                    array::xpay(Au[i], lap_u, lap_weight);
                }
                if (use_mag) { array::xpay(Au[i], u[i], rho); }
            }
            return Au;
        };

        for (size_t iter = 0; iter < outer_max; ++iter) {

            {
                // rhs := R^T y + rho [mag] (z - b) + mu [tv] ∇^T(d - bt)
                auto rhs = yT.clone();
                if (use_mag) {
                    vec_xpay(rhs, z, rho);
                    vec_xpay(rhs, b, -rho);
                }
                if (use_tv) {
                    for (size_t i = 0; i < 3; ++i) {
                        VecArray<T> d_b = d[i] - bt[i];
                        auto neg_div_db = neg_divergence(d_b);
                        array::xpay(rhs[i], neg_div_db, mu);
                    }
                }

                // x-update
                x = cgsolver(Ap, rhs, x_old, inner_max, tol, tol,
                             mag.zero_outside ? mask : no_support, (T)0, nullptr);
            }

            // z-update: project (x + b) onto {|z| = M0 inside the mask}; M0 is the
            // mean magnitude of (x + b) over the mask unless it is held fixed
            T dual_mag = 0;
            if (use_mag) {
                VecArray<T> z_new = x.clone();
                z_new += b;
                if (mag.estimate_Ms || mag.Ms <= T(0)) {
                    M0 = mean_magnitude(z_new, roi);
                }
                project_magnitude_inplace(z_new, M0, roi, &z, 0.0, mag.zero_outside);
                if (history) { dual_mag = rho * (z_new - z).norm2(); }
                z = std::move(z_new);
                b += x;
                b -= z;
            }

            // d-update: isotropic TV shrinkage of (grad x + bt)
            if (use_tv) {
                std::array<VecArray<T>, 3> dx;
                for (size_t i = 0; i < 3; ++i) { dx[i] = grad(x[i]); }
                auto sk = compute_sk(dx, bt);
                for (size_t i = 0; i < 3; ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        d[i][j] = shrink(dx[i][j], bt[i][j], sk, w.lambda_tv, mu);
                        bt[i][j] += dx[i][j] - d[i][j];
                    }
                }
            }

            // convergence: relative change of x (and ||x - z|| if const_mag)
            T xnorm = x.norm2() + (T)MAG_EPSILON;
            T rel_dx = (x - x_old).norm2() / (x_old.norm2() + (T)MAG_EPSILON);
            bool converged = rel_dx < xtol;
            T rel_r = 0;
            if (use_mag) {
                rel_r = (x - z).norm2() / xnorm;
                converged = converged && rel_r < xtol;
            }
            if (history) {
                AdmmIter<T> h;
                h.iter = iter;
                h.primal_mag = rel_r;
                h.dual_mag = dual_mag / xnorm;
                h.rho = rho;
                h.Ms = M0;
                history->iters.push_back(h);
            }
            if (logger) {
                if (use_mag) {
                    logger->log(
                        std::format("Outer iter: {}, Ms: {:.6e}, rho: {:.6e}, "
                                    "||dx||/||x||: {:.6e}, ||x-z||/||x||: {:.6e}\n",
                                    iter, M0, rho, rel_dx, rel_r));
                } else {
                    logger->log(std::format("Outer iter: {}, ||dx||/||x||: {:.6e}\n",
                                            iter, rel_dx));
                }
            }
            if (history && use_mag && (converged || iter + 1 == outer_max)) {
                history->m = x.clone();
            }
            x_old = std::move(x);
            if (converged) { break; }

            // rho continuation (monotone). The unscaled multiplier rho * b is
            // kept fixed by rescaling b.
            if (use_mag && mag.rho_growth > T(1)) {
                T rho_new = rho * mag.rho_growth;
                if (mag.rho_max > T(0)) { rho_new = std::min(rho_new, mag.rho_max); }
                if (rho_new != rho) {
                    b *= rho / rho_new;
                    rho = rho_new;
                }
            }
        }
        if (history) {
            history->rho = rho;
            history->Ms = M0;
            if (use_mag) { history->b_z = b.clone(); }
        }
        // x_old holds the last x; return the feasible variable if constrained
        return use_mag ? std::move(z) : std::move(x_old);
    }

    // explicit template instantiations for float and double
    template float mean_magnitude(const VecArray<float> &);
    template double mean_magnitude(const VecArray<double> &);
    template VecArray<float> project_magnitude(const VecArray<float> &, float);
    template VecArray<double> project_magnitude(const VecArray<double> &, double);
    template VecArray<float>
    constmag_l2(const gpuFunction<float> &A, const VecArray<float> &yT,
                const VecArray<float> &x0, float lambda_l2, float mu,
                size_t outer_max, size_t inner_max, float tol, float xtol,
                const DeviceArray<float> &mask, Logger *logger, bool zero_outside);
    template VecArray<double>
    constmag_l2(const gpuFunction<double> &A, const VecArray<double> &yT,
                const VecArray<double> &x0, double lambda_l2, double mu,
                size_t outer_max, size_t inner_max, double tol, double xtol,
                const DeviceArray<double> &mask, Logger *logger, bool zero_outside);
    template VecArray<float>
    regularized_admm(const gpuFunction<float> &, const VecArray<float> &,
                     const VecArray<float> &, const RegWeights<float> &, float,
                     size_t, size_t, float, float, const DeviceArray<float> &,
                     Logger *, const MagOptions<float> &, AdmmHistory<float> *);
    template VecArray<double>
    regularized_admm(const gpuFunction<double> &, const VecArray<double> &,
                     const VecArray<double> &, const RegWeights<double> &, double,
                     size_t, size_t, double, double, const DeviceArray<double> &,
                     Logger *, const MagOptions<double> &, AdmmHistory<double> *);
    template VecArray<float> project_sphere(const VecArray<float> &, float,
                                            const VecArray<float> &, float);
    template VecArray<double> project_sphere(const VecArray<double> &, double,
                                             const VecArray<double> &, double);
} // namespace tomocam::gpu::opt
