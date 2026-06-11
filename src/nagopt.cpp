/* -------------------------------------------------------------------------------
 * Tomocam Copyright (c) 2018
 *
 * The Regents of the University of California, through Lawrence Berkeley
 *National Laboratory (subject to receipt of any required approvals from the
 *U.S. Dept. of Energy). All rights reserved.
 *
 * If you have questions about your rights to use or distribute this software,
 * please contact Berkeley Lab's Innovation & Partnerships Office at
 *IPO@lbl.gov.
 *
 * NOTICE. This Software was developed under funding from the U.S. Department of
 * Energy and the U.S. Government consequently retains certain rights. As such,
 *the U.S. Government has been granted for itself and others acting on its
 *behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software
 *to reproduce, distribute copies to the public, prepare derivative works, and
 * perform publicly and display publicly, and to permit other to do so.
 *---------------------------------------------------------------------------------
 */

#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <tuple>

#include "array_ops.h"
#include "logger.h"
#include "optimize.h"

namespace tomocam::opt {
    template <typename T>
    Array<T> nagopt(const std::function<Array<T>(const Array<T> &)> &grad,
                    const std::function<T(const Array<T> &)> &loss, Array<T> &x,
                    size_t max_iters, T lipschitz, T tol, T xtol,
                    size_t max_inner_iters, Logger *logger) {

        Logger default_logger(LogMode::STDOUT);
        if (!logger) logger = &default_logger;

        logger->log(
            std::format("NAG Optimization started: max_iters={}, lipschitz={}, "
                        "tol={:.2e}, xtol={:.2e}\n",
                        max_iters, lipschitz, tol, xtol));

        // initialize
        Array<T> xold = x.clone();

        T t = 1;
        T tnew = 1;

        // gradient step size
        T step = 0.1 / lipschitz;
        T step_prev = step;
        T xerr = std::numeric_limits<T>::max();

        for (size_t iter = 0; iter < max_iters; iter++) {

            // line search
            size_t inner_iter_final = 0;
            for (size_t inner_iter = 0; inner_iter < max_inner_iters; inner_iter++) {
                inner_iter_final = inner_iter;

                // update theta
                T beta = tnew * (1 / t - 1);
                tnew = 0.5 * (std::sqrt(std::pow(t, 4) + 4 * std::pow(t, 2)) -
                              std::pow(t, 2));

                // update y
                auto y = x + (x - xold) * beta;
                auto g = grad(y);

                // update x
                x = y - g * step;

                // check if step size is small enough
                T ex = loss(x);
                T ey = loss(y);
                T gy = 0.5 * step * array::norm2(g);
                if (ex > (ey + gy))
                    step *= 0.9;
                else {
                    step = step_prev;
                    t = tnew;
                    xerr = array::norm2(x - xold) / array::norm2(x);
                    if (xerr < xtol) {
                        logger->log(std::format(
                            "Convergence achieved at iter {}: x-error {:.6e} < xtol "
                            "{:.2e}\n",
                            iter, xerr, xtol));
                        return x.clone();
                    }
                    xold = x.clone();
                    break;
                }
            }
            T e = loss(x);
            logger->log(
                std::format("Iter {:3d}: loss={:.6e}, x-error={:.6e}, step={:.6e}, "
                            "inner_iters={}\n",
                            iter, e, xerr, step, inner_iter_final + 1));
            // check convergence
            if (e < tol) {
                logger->log(std::format(
                    "Convergence achieved at iter {}: loss {:.6e} < tol {:.2e}\n",
                    iter, e, tol));
                break;
            }
        }
        return x.clone();
    }

    template <typename T>
    VecArray<T> nagopt(const std::function<VecArray<T>(const VecArray<T> &)> &grad,
                       const std::function<T(const VecArray<T> &)> &loss,
                       VecArray<T> &x, size_t max_iters, T lipschitz, T tol, T xtol,
                       size_t max_inner_iters, Logger *logger) {

        Logger default_logger(LogMode::STDOUT);
        if (!logger) logger = &default_logger;

        logger->log(std::format("NAG Optimization (vector) started: max_iters={}, "
                                "lipschitz={}, tol={:.2e}, xtol={:.2e}\n",
                                max_iters, lipschitz, tol, xtol));

        // initialize
        VecArray<T> xold = x.clone();

        T t = 1;
        T tnew = 1;

        // gradient step size
        T step = 0.1 / lipschitz;
        T step_prev = step;
        T xerr = std::numeric_limits<T>::max();

        for (size_t iter = 0; iter < max_iters; iter++) {

            // line search
            size_t inner_iter_final = 0;
            for (size_t inner_iter = 0; inner_iter < max_inner_iters; inner_iter++) {
                inner_iter_final = inner_iter;

                // update theta
                T beta = tnew * (1 / t - 1);
                tnew = 0.5 * (std::sqrt(std::pow(t, 4) + 4 * std::pow(t, 2)) -
                              std::pow(t, 2));

                // update y and x
                auto y = x + (x - xold) * beta;
                auto g = grad(y);
                x = y - g * step;

                // check if step size is small enough
                T ex = loss(x);
                T ey = loss(y);
                T gy = 0.5 * step * g.norm2();
                if (ex > (ey + gy))
                    step *= 0.9;
                else {
                    step = step_prev;
                    t = tnew;
                    xerr = (x - xold).norm2() / x.norm2();
                    if (xerr < xtol) {
                        logger->log(std::format(
                            "Convergence achieved at iter {}: x-error {:.6e} < xtol "
                            "{:.2e}\n",
                            iter, xerr, xtol));
                        return x.clone();
                    }
                    xold = x.clone();
                    break;
                }
            }
            T e = loss(x);
            logger->log(std::format("Iter {:3d}: loss={:.6e}, x-error={:.6e}\n",
                                    iter, e, xerr));
            // check convergence
            if (e < tol) {
                logger->log(std::format(
                    "Convergence achieved at iter {}: loss {:.6e} < tol {:.2e}\n",
                    iter, e, tol));
                break;
            }
        }
        return x.clone();
    }

    // explicit template instantiation
    template Array<float>
    nagopt<float>(const std::function<Array<float>(const Array<float> &)> &,
                  const std::function<float(const Array<float> &)> &, Array<float> &,
                  size_t, float, float, float, size_t, Logger *);
    template Array<double>
    nagopt<double>(const std::function<Array<double>(const Array<double> &)> &,
                   const std::function<double(const Array<double> &)> &,
                   Array<double> &, size_t, double, double, double, size_t,
                   Logger *);

    template VecArray<float>
    nagopt<float>(const std::function<VecArray<float>(const VecArray<float> &)> &,
                  const std::function<float(const VecArray<float> &)> &,
                  VecArray<float> &, size_t, float, float, float, size_t, Logger *);
    template VecArray<double>
    nagopt<double>(const std::function<VecArray<double>(const VecArray<double> &)> &,
                   const std::function<double(const VecArray<double> &)> &,
                   VecArray<double> &, size_t, double, double, double, size_t,
                   Logger *);

} // namespace tomocam::opt
