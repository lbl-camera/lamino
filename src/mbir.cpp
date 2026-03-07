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
#include <cassert>
#include <format>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "optimize.h"
#include "padding.h"
#include "polar_grid.h"
#include "tomocam.h"

namespace tomocam {

    template <typename T>
    std::array<Array<T>, 3> MBIR(const Array<T> &proj, const std::vector<T> &angles,
                                 T gamma, const dims_t &recon_dims, size_t max_iter,
                                 T sigma, T p, T tol, T xtol) {

        // padding factor
        constexpr double PAD_FACTOR = 1.42;
        T padding = static_cast<T>(PAD_FACTOR);

        // adjust reconstruction dimensions
        dims_t out_dims = recon_dims;
        out_dims.n1 = static_cast<size_t>(recon_dims.n1 * padding);
        if (out_dims.n1 % 2 == 0) {
            out_dims.n1 -= 1; // make sure n1 is odd
        }

        out_dims.n2 = static_cast<size_t>(recon_dims.n2 * padding);
        if (out_dims.n2 % 2 == 0) {
            out_dims.n2 -= 1; // make sure n2 is odd
        }
        out_dims.n3 = static_cast<size_t>(recon_dims.n3 * padding);
        if (out_dims.n3 % 2 == 0) {
            out_dims.n3 -= 1; // make sure n3 is odd
        }

        // normalize projections
        T proj_max = array::max(proj);
        auto y = proj / proj_max;

        // zero-pad projections by sqrt(2) to avoid aliasing
        y = pad2d(y, padding, PadType::SYMMETRIC);

        // setup polar grid
        size_t nrows = y.nrows();
        size_t ncols = y.ncols();
        auto polar_grid = PolarGrid<T>(angles, nrows, ncols, gamma);

        // backproject measurements to get yT
        auto yT = adjoint(y, polar_grid, out_dims, gamma);

        // setup gradient operator
        std::function<std::array<Array<T>, 3>(const std::array<Array<T>, 3> &)>
            grad = [&](const std::array<Array<T>, 3> &x) {
                // gradient data
                auto g_data = gradient(x, yT, polar_grid, gamma);

                // apply qggmrf penalty
                for (size_t i = 0; i < 3; ++i) {
                    opt::qggmrf(x[i], g_data[i], sigma, p);
                }
                return g_data;
            };

        // setup loss function
        std::function<T(const std::array<Array<T>, 3> &)> loss =
            [&](const std::array<Array<T>, 3> &x) {
                return array::norm2(forward(x, polar_grid, gamma) - y);
            };

        // lipshitz constant
        std::array<Array<T>, 3> xtmp;
        for (size_t i = 0; i < 3; ++i) { xtmp[i] = Array<T>::ones(out_dims); }
        auto Axtmp = forward(xtmp, polar_grid, gamma);
        auto gtmp = adjoint(Axtmp, polar_grid, out_dims, gamma);
        T L = 0;
        for (size_t i = 0; i < 3; ++i) {
            // estimate lipshitz constant
            L = std::max(L, array::max(array::abs(gtmp[i])));
        }

        // initial guess
        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = Array<T>::ones(out_dims) * 0.9; }

        auto recon_m = opt::nagopt(grad, loss, x0, max_iter, L, tol, xtol);

        // crop to original dimensions
        std::array<Array<T>, 3> recon_magnetisation;
        for (size_t i = 0; i < 3; ++i) {
            recon_magnetisation[i] =
                crop3d(recon_m[i], recon_dims, PadType::SYMMETRIC);
        }
        return recon_magnetisation;
    }

    // Explicit template instantiations
    template std::array<Array<float>, 3> MBIR(const

} // namespace tomocam
