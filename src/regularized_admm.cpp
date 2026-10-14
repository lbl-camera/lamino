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
#include <cmath>
#include <execution>
#include <format>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "bregman.h"
#include "optimize.h"

namespace tomocam::opt {

    constexpr double MAG_EPSILON = 1e-12;
    constexpr double TV_EPSILON = 1e-8; // same as split_bregman

    namespace {
        // Number of voxels in the support mask (mask > 0.5). The mask must have
        // the volume dims; voxels outside it are padding and are left
        // unconstrained.
        template <typename T>
        size_t support_size(const Array<T> &mask, const dims_t &dims) {
            if (mask.dims() != dims) {
                throw std::invalid_argument(
                    "constant-magnitude mask must have the volume dimensions");
            }
            size_t n =
                static_cast<size_t>(std::count_if(mask.begin(), mask.end(), [](T v) {
                    return v > static_cast<T>(0.5);
                }));
            if (n == 0) {
                throw std::invalid_argument("constant-magnitude mask is empty");
            }
            return n;
        }

        // copy of x0, zeroed outside the support if zero_outside
        template <typename T>
        VecArray<T> initial_guess(const VecArray<T> &x0, const Array<T> &mask,
                                  bool zero_outside) {
            VecArray<T> x = clone(x0);
            if (zero_outside) {
                for (size_t i = 0; i < 3; ++i) { x[i] *= mask; }
            }
            return x;
        }

        template <typename T>
        T mean_magnitude(const VecArray<T> &v, const Array<T> &mask, size_t count) {
            size_t n = v[0].size();
            const T *px = v[0].begin(), *py = v[1].begin(), *pz = v[2].begin();
            const T *pm = mask.begin();
            std::vector<size_t> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            double sum = std::transform_reduce(
                std::execution::par_unseq, idx.begin(), idx.end(), 0.0,
                std::plus<double>(), [=](size_t i) {
                    if (!(pm[i] > T(0.5))) return 0.0;
                    double x = px[i], y = py[i], z = pz[i];
                    return std::sqrt(x * x + y * y + z * z);
                });
            return static_cast<T>(sum / static_cast<double>(count));
        }

        // Project v onto |z| = M0 inside the mask. A voxel with |v| <= eps keeps
        // its previous value (the projection is not unique there); with no
        // z_prev it becomes zero. Outside the mask z = v, or 0 if zero_outside.
        template <typename T>
        VecArray<T>
        project_magnitude(const VecArray<T> &v, T M0, const Array<T> &mask,
                          const VecArray<T> *z_prev = nullptr,
                          bool zero_outside = false, double eps = MAG_EPSILON) {
            VecArray<T> rv = clone(v);
            size_t n = rv[0].size();
            T *px = rv[0].begin(), *py = rv[1].begin(), *pz = rv[2].begin();
            const T *pm = mask.begin();
            const T *qx = z_prev ? (*z_prev)[0].begin() : nullptr;
            const T *qy = z_prev ? (*z_prev)[1].begin() : nullptr;
            const T *qz = z_prev ? (*z_prev)[2].begin() : nullptr;
            std::vector<size_t> idx(n);
            std::iota(idx.begin(), idx.end(), 0);
            std::for_each(std::execution::par_unseq, idx.begin(), idx.end(),
                          [=](size_t i) {
                              if (!(pm[i] > T(0.5))) {
                                  if (zero_outside) { px[i] = py[i] = pz[i] = T(0); }
                                  return;
                              }
                              double x = px[i], y = py[i], z = pz[i];
                              double mag = std::sqrt(x * x + y * y + z * z);
                              if (mag > eps) {
                                  double s = M0 / mag;
                                  px[i] = static_cast<T>(x * s);
                                  py[i] = static_cast<T>(y * s);
                                  pz[i] = static_cast<T>(z * s);
                              } else if (qx) {
                                  px[i] = qx[i];
                                  py[i] = qy[i];
                                  pz[i] = qz[i];
                              } else {
                                  px[i] = py[i] = pz[i] = T(0);
                              }
                          });
            return rv;
        }
    } // namespace

    template <typename T>
    T mean_magnitude(const VecArray<T> &v) {
        auto mask = Array<T>::ones(v[0].dims());
        return mean_magnitude(v, mask, mask.size());
    }

    template <typename T>
    VecArray<T> project_magnitude(const VecArray<T> &v, T M0) {
        return project_magnitude(v, M0, Array<T>::ones(v[0].dims()));
    }

    template <typename T>
    VecArray<T> project_sphere(const VecArray<T> &v, T Ms, const VecArray<T> &z_prev,
                               T eps) {
        return project_magnitude(v, Ms, Array<T>::ones(v[0].dims()), &z_prev, false,
                                 static_cast<double>(eps));
    }

    template <typename T>
    VecArray<T>
    constmag_l2(const Function<T> &A, const VecArray<T> &yT, const VecArray<T> &x0,
                T lambda_l2, T mu, size_t outer_max, size_t inner_max, T tol, T xtol,
                const Array<T> &mask, Logger *logger, bool zero_outside) {

        if (x0[0].dims() != yT[0].dims()) {
            throw std::invalid_argument(
                "yT is not the raw data, but backprojected data, so it should have "
                "the same dimensions as x0");
        }
        auto dims = x0[0].dims();
        const size_t nmask = support_size(mask, dims);
        const Array<T> no_support;

        // the magnitude is uniform but unknown: start from the mean magnitude
        // of the initial guess and re-estimate it in every z-update
        T M0 = mean_magnitude(x0, mask, nmask);
        if (logger)
            logger->log(std::format("Initial constant magnitude M0: {:.6e}\n", M0));

        VecArray<T> x = initial_guess(x0, mask, zero_outside);
        VecArray<T> x_old = clone(x);
        VecArray<T> z = project_magnitude(
            x, M0, mask, static_cast<const VecArray<T> *>(nullptr), zero_outside);
        VecArray<T> b = {Array<T>::zeros(dims), Array<T>::zeros(dims),
                         Array<T>::zeros(dims)};

        // Ap = (A^TA + lambda_l2 ∇^T∇ + mu I) u
        Function<T> Ap = [&](const VecArray<T> &u) {
            auto du = A(u);
            for (size_t i = 0; i < 3; ++i) {
                auto lap = laplacian<T>(u[i]);
                lap *= lambda_l2;
                du[i] += lap;
                auto ui = u[i].clone();
                ui *= mu;
                du[i] += ui;
            }
            return du;
        };

        for (size_t iter = 0; iter < outer_max; ++iter) {

            // rhs := R^T y + mu (z - b)
            VecArray<T> rhs;
            for (size_t i = 0; i < 3; ++i) {
                auto t = z[i] - b[i];
                t *= mu;
                rhs[i] = yT[i] + t;
            }

            // x-update
            x = cgsolver<T>(Ap, rhs, x_old, inner_max, tol, tol,
                            zero_outside ? mask : no_support, T(0), nullptr);

            // z-update: project (x + b) onto {|z| = M0 inside the mask, M0 free}.
            // The optimal M0 is the mean magnitude of (x + b) over the mask;
            // outside the mask z = x + b, so b stays zero there and the padding
            // is left unconstrained.
            auto v = x + b;
            M0 = mean_magnitude(v, mask, nmask);
            z = project_magnitude(v, M0, mask,
                                  static_cast<const VecArray<T> *>(nullptr),
                                  zero_outside);

            // dual update
            for (size_t i = 0; i < 3; ++i) {
                b[i] += x[i];
                b[i] -= z[i];
            }

            // convergence: relative change of x and primal residual ||x - z||
            T xnorm = std::sqrt(dot(x, x)) + static_cast<T>(MAG_EPSILON);
            auto dx = x - x_old;
            auto r = x - z;
            T rel_dx = std::sqrt(dot(dx, dx)) /
                       (std::sqrt(dot(x_old, x_old)) + static_cast<T>(MAG_EPSILON));
            T rel_r = std::sqrt(dot(r, r)) / xnorm;
            if (logger)
                logger->log(std::format("Outer iter: {}, M0: {:.6e}, "
                                        "||dx||/||x||: {:.6e}, "
                                        "||x-z||/||x||: {:.6e}\n",
                                        iter, M0, rel_dx, rel_r));
            x_old = clone(x);
            if (rel_dx < xtol && rel_r < xtol) { break; }
        }
        // return the feasible variable
        return z;
    }

    template <typename T>
    VecArray<T> regularized_admm(const Function<T> &A, const VecArray<T> &yT,
                                 const VecArray<T> &x0, const RegWeights<T> &w, T mu,
                                 size_t outer_max, size_t inner_max, T tol, T xtol,
                                 const Array<T> &mask, Logger *logger,
                                 const MagOptions<T> &mag, AdmmHistory<T> *history) {

        if (x0[0].dims() != yT[0].dims()) {
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
        auto dims = x0[0].dims();
        const size_t nmask = support_size(mask, dims);
        const Array<T> no_support;
        const T eps = static_cast<T>(MAG_EPSILON);

        VecArray<T> x = initial_guess(x0, mask, mag.zero_outside);
        VecArray<T> x_old = clone(x);

        // magnitude splitting: z = x, dual b, penalty rho (continued upwards)
        T rho = mag.rho > T(0) ? mag.rho : mu;
        T M0 = 0;
        VecArray<T> z, z_old, b;
        if (use_mag) {
            M0 = (!mag.estimate_Ms && mag.Ms > T(0))
                     ? mag.Ms
                     : mean_magnitude(x0, mask, nmask);
            if (logger)
                logger->log(std::format("Initial constant magnitude Ms: {:.6e}, "
                                        "rho: {:.6e}\n",
                                        M0, rho));
            // z0 = P(x0); zero vectors in the mask get the default direction +x
            VecArray<T> dflt = {Array<T>::ones(dims) * M0, Array<T>::zeros(dims),
                                Array<T>::zeros(dims)};
            z = project_magnitude(x, M0, mask, &dflt, mag.zero_outside);
            b = {Array<T>::zeros(dims), Array<T>::zeros(dims),
                 Array<T>::zeros(dims)};
            if (mag.zero_outside) {
                // z0 is 0 outside the mask, so start with the matching multiplier
                for (size_t i = 0; i < 3; ++i) { b[i] = x[i] - z[i]; }
            }
        }

        // TV splitting: d = grad x, dual bt (3x3 matrices)
        std::array<std::array<Array<T>, 3>, 3> d, bt;
        if (use_tv) {
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    d[i][j] = Array<T>::zeros(dims);
                    bt[i][j] = Array<T>::zeros(dims);
                }
            }
        }

        // Ap = (A^TA + (lambda_l2 + mu [tv]) ∇^T∇ + rho [mag] I) u
        // beta only enters the system matrix, never the right-hand side. rho is
        // captured by reference: it may grow between outer iterations.
        const T lap_weight = w.lambda_l2 + (use_tv ? mu : T(0));
        const bool use_lap = use_tv || w.l2();
        Function<T> Ap = [&](const VecArray<T> &u) {
            auto du = A(u);
            for (size_t i = 0; i < 3; ++i) {
                if (use_lap) {
                    auto lap = laplacian<T>(u[i]);
                    lap *= lap_weight;
                    du[i] += lap;
                }
                if (use_mag) {
                    auto ui = u[i].clone();
                    ui *= rho;
                    du[i] += ui;
                }
            }
            return du;
        };

        auto sumsq = [](const Array<T> &a) {
            return std::transform_reduce(a.begin(), a.end(), T(0), std::plus<T>(),
                                         [](T v) { return v * v; });
        };
        auto norm = [](const VecArray<T> &v) { return std::sqrt(dot(v, v)); };

        for (size_t iter = 0; iter < outer_max; ++iter) {

            // rhs := R^T y + rho [mag] (z - b) + mu [tv] ∇^T(d - bt)
            VecArray<T> rhs = clone(yT);
            for (size_t i = 0; i < 3; ++i) {
                if (use_mag) {
                    auto t = z[i] - b[i];
                    t *= rho;
                    rhs[i] += t;
                }
                if (use_tv) { rhs[i] += divergence(d[i] - bt[i]) * mu; }
            }

            // x-update
            x = cgsolver<T>(Ap, rhs, x_old, inner_max, tol, tol,
                            mag.zero_outside ? mask : no_support, T(0), nullptr);

            // z-update: project (x + b) onto {|z| = M0 inside the mask}. M0 is the
            // mean magnitude of (x + b) over the mask unless it is held fixed.
            // Outside the mask z = x + b (b stays zero there: padding is free)
            // or z = 0 if zero_outside.
            T primal_mag = 0, dual_mag = 0;
            if (use_mag) {
                z_old = history ? clone(z) : VecArray<T>();
                auto v = x + b;
                if (mag.estimate_Ms || mag.Ms <= T(0)) {
                    M0 = mean_magnitude(v, mask, nmask);
                }
                z = project_magnitude(v, M0, mask, &z, mag.zero_outside);
                for (size_t i = 0; i < 3; ++i) {
                    b[i] += x[i];
                    b[i] -= z[i];
                }
            }

            // d-update: isotropic TV shrinkage of (grad x + bt)
            T primal_tv = 0, dual_tv = 0;
            std::array<std::array<Array<T>, 3>, 3> grad_x;
            if (use_tv) {
                std::array<std::array<Array<T>, 3>, 3> d_old;
                if (history) {
                    for (size_t i = 0; i < 3; ++i) {
                        for (size_t j = 0; j < 3; ++j) {
                            d_old[i][j] = d[i][j].clone();
                        }
                    }
                }
                for (size_t i = 0; i < 3; ++i) { grad_x[i] = grad_u(x[i]); }

                auto sk = Array<T>::ones(dims) * static_cast<T>(TV_EPSILON);
                for (size_t i = 0; i < 3; ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        sk += (grad_x[i][j] + bt[i][j]) * (grad_x[i][j] + bt[i][j]);
                    }
                }
                std::transform(std::execution::par_unseq, sk.begin(), sk.end(),
                               sk.begin(), [](T val) { return std::sqrt(val); });

                const T lambda = w.lambda_tv;
                for (size_t i = 0; i < 3; ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        auto temp = (grad_x[i][j] + bt[i][j]) / sk;
                        std::transform(
                            std::execution::par_unseq, temp.begin(), temp.end(),
                            sk.begin(), d[i][j].begin(),
                            [lambda, mu](T val, T sk_val) {
                                return std::max(sk_val - lambda / mu, T(0)) * val;
                            });
                        bt[i][j] += grad_x[i][j] - d[i][j];
                    }
                }

                if (history) {
                    // primal ||grad x - d||, dual mu ||grad^T (d_new - d_old)||
                    T p2 = 0, d2 = 0;
                    for (size_t i = 0; i < 3; ++i) {
                        for (size_t j = 0; j < 3; ++j) {
                            auto r = grad_x[i][j] - d[i][j];
                            p2 += sumsq(r);
                        }
                        auto q = divergence(d[i] - d_old[i]);
                        d2 += sumsq(q);
                    }
                    primal_tv = std::sqrt(p2);
                    dual_tv = mu * std::sqrt(d2);
                }
            }

            // convergence: relative change of x (and ||x - z|| if const_mag)
            auto dx = x - x_old;
            T xnorm = norm(x) + eps;
            T rel_dx = norm(dx) / (norm(x_old) + eps);
            bool converged = rel_dx < xtol;
            T rel_r = 0;
            if (use_mag) {
                auto r = x - z;
                rel_r = norm(r) / xnorm;
                converged = converged && rel_r < xtol;
                if (history) {
                    primal_mag = rel_r;
                    dual_mag = rho * norm(z - z_old) / xnorm;
                }
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

            if (history) {
                AdmmIter<T> h;
                h.iter = iter;
                h.primal_tv = primal_tv / xnorm;
                h.dual_tv = dual_tv / xnorm;
                h.primal_mag = primal_mag;
                h.dual_mag = dual_mag;
                h.rho = rho;
                h.Ms = M0;
                // objective terms at x: 0.5<x,A^TAx> - <x,A^Ty>, TV,
                // lambda_l2/2|grad x|^2
                auto Ax = A(x);
                h.data_shifted = T(0.5) * dot(x, Ax) - dot(x, yT);
                T tv = 0, g2 = 0;
                if (!use_tv) {
                    for (size_t i = 0; i < 3; ++i) { grad_x[i] = grad_u(x[i]); }
                }
                auto s = Array<T>::zeros(dims);
                for (size_t i = 0; i < 3; ++i) {
                    for (size_t j = 0; j < 3; ++j) {
                        s += grad_x[i][j] * grad_x[i][j];
                    }
                }
                for (auto it = s.begin(); it != s.end(); ++it) {
                    g2 += *it;
                    tv += std::sqrt(*it);
                }
                h.tv = w.lambda_tv * tv;
                h.l2 = T(0.5) * w.lambda_l2 * g2;
                history->iters.push_back(h);
            }

            x_old = clone(x);
            if (converged) { break; }

            // rho continuation (monotone). The unscaled multiplier rho * b is
            // kept fixed by rescaling b.
            if (use_mag && mag.rho_growth > T(1)) {
                T rho_new = rho * mag.rho_growth;
                if (mag.rho_max > T(0)) { rho_new = std::min(rho_new, mag.rho_max); }
                if (rho_new != rho) {
                    for (size_t i = 0; i < 3; ++i) { b[i] *= rho / rho_new; }
                    rho = rho_new;
                }
            }
        }

        if (history) {
            history->rho = rho;
            history->Ms = M0;
            if (use_mag) {
                history->m = clone(x);
                history->b_z = clone(b);
            }
        }
        // return the feasible variable if the constraint is active
        if (use_mag) { return z; }
        return x;
    }

    template VecArray<float>
    regularized_admm(const Function<float> &, const VecArray<float> &,
                     const VecArray<float> &, const RegWeights<float> &, float,
                     size_t, size_t, float, float, const Array<float> &, Logger *,
                     const MagOptions<float> &, AdmmHistory<float> *);
    template VecArray<double>
    regularized_admm(const Function<double> &, const VecArray<double> &,
                     const VecArray<double> &, const RegWeights<double> &, double,
                     size_t, size_t, double, double, const Array<double> &, Logger *,
                     const MagOptions<double> &, AdmmHistory<double> *);
    template VecArray<float> project_sphere(const VecArray<float> &, float,
                                            const VecArray<float> &, float);
    template VecArray<double> project_sphere(const VecArray<double> &, double,
                                             const VecArray<double> &, double);

    template float mean_magnitude(const VecArray<float> &);
    template double mean_magnitude(const VecArray<double> &);
    template VecArray<float> project_magnitude(const VecArray<float> &, float);
    template VecArray<double> project_magnitude(const VecArray<double> &, double);
    template VecArray<float>
    constmag_l2(const Function<float> &A, const VecArray<float> &yT,
                const VecArray<float> &x0, float lambda_l2, float mu,
                size_t outer_max, size_t inner_max, float tol, float xtol,
                const Array<float> &mask, Logger *logger, bool zero_outside);
    template VecArray<double>
    constmag_l2(const Function<double> &A, const VecArray<double> &yT,
                const VecArray<double> &x0, double lambda_l2, double mu,
                size_t outer_max, size_t inner_max, double tol, double xtol,
                const Array<double> &mask, Logger *logger, bool zero_outside);

} // namespace tomocam::opt
