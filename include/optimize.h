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

#ifndef OPTIMIZE__H
#define OPTIMIZE__H

#include <cmath>
#include <iostream>

#include "array.h"
#include "tomocam.h"

namespace tomocam::opt {

    template <typename T>
    using Function = std::function<Array<T>(const Array<T> &)>;

    template <typename T>
    using Residual = std::function<T(const Array<T> &)>;

    // bregman.cpp
    template <typename T>
    Array<T> split_bregman(const Function<T> &A, const Array<T> &y, Array<T> x0,
                           T lambda, T mu, size_t outer_max, size_t inner_max, T tol,
                           T xtol);

    // conjgrad.cpp
    template <typename T>
    Array<T> cgsolver(const Function<T> &A, const Array<T> &y, size_t max_iter,
                      T tol);

    // nagopt.cpp
    template <typename T>
    Array<T> nagopt(const Function<T> &grad, const Residual<T> &loss, Array<T> &x,
                    size_t max_iters, T lipschitz, T tol, T xtol,
                    size_t max_inner_iters = 20);

    // lipschitz.cpp
    template <typename T>
    T lipschitz(const Function<T> &grad, const Array<T> &x0, size_t max_iters = 100,
                T tol = 1e-5);

    // qggmrf.cpp
    template <typename T>
    void qggmrf(const Array<T> &x, Array<T> &g, T sigma, T p);

} // namespace tomocam::opt

#endif // TOMOCAM_OPTIMIZE__H
