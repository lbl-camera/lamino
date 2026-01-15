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
    Array<T> nagopt(const Function<T> &grad, const Residual<T> &loss, Array<T> &x,
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
    std::array<Array<T>, 3> nagopt(
        const std::function<std::array<Array<T>, 3>(const std::array<Array<T>, 3> &)>
            &grad,
        const std::function<T(const std::array<Array<T>, 3> &)> &loss,
        std::array<Array<T>, 3> &x, size_t max_iters, T lipschitz, T tol, T xtol,
        size_t max_inner_iters, Logger *logger) {

        Logger default_logger(LogMode::STDOUT);
        if (!logger) logger = &default_logger;

        logger->log(std::format("NAG Optimization (vector) started: max_iters={}, "
                                "lipschitz={}, tol={:.2e}, xtol={:.2e}\n",
                                max_iters, lipschitz, tol, xtol));

        // initialize
        std::array<Array<T>, 3> xold;
        for (size_t i = 0; i < 3; ++i) { xold[i] = x[i].clone(); }

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
                std::array<Array<T>, 3> y;
                for (size_t i = 0; i < 3; ++i) {
                    y[i] = x[i] + (x[i] - xold[i]) * beta;
                }
                auto g = grad(y);

                // update x
                for (size_t i = 0; i < 3; ++i) { x[i] = y[i] - g[i] * step; }

                // check if step size is small enough
                T ex = loss(x);
                T ey = loss(y);
                T gy = 0;
                for (size_t i = 0; i < 3; ++i) { gy += array::norm2(g[i]); }
                gy *= 0.5 * step;
                if (ex > (ey + gy))
                    step *= 0.9;
                else {
                    step = step_prev;
                    t = tnew;
                    T xerr_sum = 0;
                    T xnorm_sum = 0;
                    for (size_t i = 0; i < 3; ++i) {
                        xerr_sum += array::norm2(x[i] - xold[i]);
                        xnorm_sum += array::norm2(x[i]);
                    }
                    xerr = xerr_sum / xnorm_sum;
                    if (xerr < xtol) {
                        logger->log(std::format(
                            "Convergence achieved at iter {}: x-error {:.6e} < xtol "
                            "{:.2e}\n",
                            iter, xerr, xtol));
                        std::array<Array<T>, 3> result;
                        for (size_t i = 0; i < 3; ++i) { result[i] = x[i].clone(); }
                        return result;
                    }
                    for (size_t i = 0; i < 3; ++i) { xold[i] = x[i].clone(); }
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
        std::array<Array<T>, 3> result;
        for (size_t i = 0; i < 3; ++i) { result[i] = x[i].clone(); }
        return result;
    }

    // explicit template instantiation
    template Array<float> nagopt<float>(const Function<float> &,
                                        const Residual<float> &, Array<float> &,
                                        size_t, float, float, float, size_t,
                                        Logger *);
    template Array<double> nagopt<double>(const Function<double> &,
                                          const Residual<double> &, Array<double> &,
                                          size_t, double, double, double, size_t,
                                          Logger *);

    template std::array<Array<float>, 3>
    nagopt<float>(const std::function<std::array<Array<float>, 3>(
                      const std::array<Array<float>, 3> &)> &,
                  const std::function<float(const std::array<Array<float>, 3> &)> &,
                  std::array<Array<float>, 3> &, size_t, float, float, float, size_t,
                  Logger *);
    template std::array<Array<double>, 3> nagopt<double>(
        const std::function<
            std::array<Array<double>, 3>(const std::array<Array<double>, 3> &)> &,
        const std::function<double(const std::array<Array<double>, 3> &)> &,
        std::array<Array<double>, 3> &, size_t, double, double, double, size_t,
        Logger *);

} // namespace tomocam::opt
