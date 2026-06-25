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
#include "projection.h"
#include "recon_params.h"

namespace tomocam {

    template <typename T>
    std::array<Array<T>, 3> MBIR(const std::vector<Dataset_t<T>> &datasets,
                                 const ReconParams &params) {

        // padding factor
        T padfac = static_cast<T>(params.PAD_FACTOR);

        // adjust reconstruction dimensions
        dims_t proj_dims = std::get<0>(datasets[0]).dims();
        dims_t output_dims = params.recon_dims;
        dims_t recon_dims = {output_dims.n1,
                             static_cast<size_t>(proj_dims.n2 * padfac),
                             static_cast<size_t>(proj_dims.n3 * padfac)};

        // array for accumulated backprojections
        std::array<Array<T>, 3> yT;
        for (size_t i = 0; i < 3; ++i) { yT[i] = Array<T>::zeros(recon_dims); }

        size_t n_datasets = datasets.size();
        std::vector<PolarGrid<T>> polar_grids(n_datasets);
        std::vector<T> gammas(n_datasets);

        T proj_max = 0.0;
        for (const auto &[proj, angles, gamma_ref] : datasets) {
            proj_max = std::max(proj_max, array::max(proj));
        }

        for (size_t j = 0; j < n_datasets; ++j) {
            auto &[proj, angles, gamma_ref] = datasets[j];
            gammas[j] = gamma_ref;

            // scale projections
            auto y = proj / proj_max;
            y = pad2d(y, padfac, PadType::SYMMETRIC);

            // setup polar grid
            size_t nrows = y.nrows();
            size_t ncols = y.ncols();
            polar_grids[j] =
                std::move(PolarGrid<T>(angles, nrows, ncols, gammas[j]));

            // backproject and accumulate yT
            auto yTmp = adjoint(y, polar_grids[j], recon_dims, gammas[j]);
            for (size_t i = 0; i < 3; ++i) { yT[i] += yTmp[i]; }
        }

        // initial guess
        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = Array<T>::zeros(recon_dims); }

        // setup linear operator
        opt::Function<T> A = [&polar_grids,
                              &gammas](const std::array<Array<T>, 3> &m) {
            std::array<Array<T>, 3> Ax = sysmat(m, polar_grids[0], gammas[0]);
            for (size_t i = 1; i < polar_grids.size(); ++i) {
                auto tmp = sysmat(m, polar_grids[i], gammas[i]);
                for (size_t j = 0; j < 3; ++j) { Ax[j] += tmp[j]; }
            }
            return Ax;
        };

        std::array<Array<T>, 3> recon_m;
        switch (params.regularizer) {
            case Regularizer::SPLIT_BREGMAN:
                recon_m = opt::split_bregman<T>(
                    A, yT, x0, params.lambda, params.mu, params.maxIters,
                    params.innerIters, params.tol, params.xtol, output_dims);
                break;
            case Regularizer::UNCONSTRAINED:
                recon_m = opt::cgsolver<T>(A, yT, x0, params.maxIters, params.tol,
                                           params.xtol, output_dims);
                // TV regularization
                break;
            default: throw std::invalid_argument("Unsupported regularizer");
        }

        // crop to original dimensions
        std::array<Array<T>, 3> recon_magnetisation;
        for (size_t i = 0; i < 3; ++i) {
            recon_magnetisation[i] =
                crop3d(recon_m[i], output_dims, PadType::SYMMETRIC);
        }
        return recon_magnetisation;
    }

    // Explicit template instantiations
    template std::array<Array<float>, 3>
    MBIR(const std::vector<Dataset_t<float>> &datasets, const ReconParams &params);
    template std::array<Array<double>, 3>
    MBIR(const std::vector<Dataset_t<double>> &datasets, const ReconParams &params);

} // namespace tomocam
