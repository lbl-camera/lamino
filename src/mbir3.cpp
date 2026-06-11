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
#include "config.h"
#include "optimize.h"
#include "padding.h"
#include "polar_grid.h"
#include "tomocam.h"

namespace tomocam {

    /** Dataset_t type definition
     * @tparam T data type
     * @brief Tuple containing projection data, projection angles, and orientation
     * gamma
     */
    template <typename T>
    using Dataset_t = std::tuple<Array<T>, std::vector<T>, T>;

    template <typename T>
    std::array<Array<T>, 3> MBIR3(const std::vector<Dataset_t<T>> &datasets,
                                  const dims_t &recon_dims,
                                  const ReconParams &recon_params) {

        T padding = static_cast<T>(recon_params.PAD_FACTOR);
        auto pad = [padding](size_t n) {
            size_t npad = 2 * (static_cast<size_t>(n * (padding - 1)) / 2);
            return n + npad;
        };

        // adjust reconstruction dimensions
        dims_t out_dims = recon_dims;
        out_dims.n1 = pad(recon_dims.n1);
        out_dims.n2 = pad(recon_dims.n2);
        out_dims.n3 = pad(recon_dims.n3);

        // setup system matrices and backprojections
        size_t n_datasets = datasets.size();
        std::array<Array<T>, 3> yT;
        for (size_t i = 0; i < 3; ++i) { yT[i] = Array<T>::zeros(out_dims); }
        std::vector<T> gamma(n_datasets);
        std::vector<PolarGrid<T>> polar_grid(n_datasets);
        for (size_t j = 0; j < n_datasets; ++j) {
            auto &[proj, angles, gamma_ref] = datasets[j];
            gamma[j] = gamma_ref;

            // normalize projections
            T proj_max = array::max(proj);
            auto y = proj / proj_max;

            // zero-pad projections by sqrt(2) to avoid aliasing
            y = pad2d(y, padding, PadType::SYMMETRIC);

            // setup polar grid
            size_t nrows = y.nrows();
            size_t ncols = y.ncols();
            polar_grid[j] = std::move(PolarGrid<T>(angles, nrows, ncols, gamma[j]));

            // backproject measurements to get yT
            auto yTmp = adjoint(y, polar_grid[j], out_dims, gamma[j]);
            for (size_t i = 0; i < 3; ++i) { yT[i] += yTmp[i]; }
        }

        // Define gradient function for NAG optimizer
        auto grad = [&polar_grid, &gamma, &yT,
                     &recon_params](const std::array<Array<T>, 3> &x) {
            auto g = gradient(x, yT, polar_grid[0], gamma[0]);
            for (size_t j = 1; j < gamma.size(); ++j) {
                auto gTmp = gradient(x, yT, polar_grid[j], gamma[j]);
                for (size_t i = 0; i < 3; ++i) { g[i] += gTmp[i]; }
            }
            // add TV regularization gradient in-place
            for (size_t i = 0; i < 3; ++i) {
                opt::qggmrf<T>(g[i], x[i], recon_params.sigma, recon_params.p);
            }
            return g;
        };

        T yTy = 0;
        for (size_t i = 0; i < 3; ++i) { yTy += array::dot(yT[i], yT[i]); }

        // Define loss function for NAG optimizer (data fidelity only)
        auto loss = [&polar_grid, &gamma, &yT, &recon_params,
                     yTy](const std::array<Array<T>, 3> &x) {
            T e = residual(x, yT, polar_grid[0], yTy, gamma[0]);
            for (size_t j = 1; j < gamma.size(); ++j) {
                e += residual(x, yT, polar_grid[j], yTy, gamma[j]);
            }
            return e;
        };

        // initial guess
        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = Array<T>::zeros(out_dims); }

        // Run NAG optimization
        auto recon_m = opt::nagopt<T>(
            grad, loss, x0, recon_params.maxIters, recon_params.mu, recon_params.tol,
            recon_params.xtol, recon_params.innerIters, nullptr);

        std::cout << "returned from optimization" << std::endl;
        // crop to original dimensions
        std::array<Array<T>, 3> recon_magnetisation;
        for (size_t i = 0; i < 3; ++i) {
            recon_magnetisation[i] =
                crop3d(recon_m[i], recon_dims, PadType::SYMMETRIC);
        }
        return recon_magnetisation;
    }

    // Explicit template instantiations
    template std::array<Array<float>, 3>
    MBIR3(const std::vector<Dataset_t<float>> &datasets, const dims_t &recon_dims,
          const ReconParams &recon_params);
    template std::array<Array<double>, 3>
    MBIR3(const std::vector<Dataset_t<double>> &datasets, const dims_t &recon_dims,
          const ReconParams &recon_params);

} // namespace tomocam
