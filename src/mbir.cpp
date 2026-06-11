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
    std::array<Array<T>, 3> MBIR(std::vector<Dataset_t<T>> datasets,
                                 ReconParams params) {

        // padding factor
        T padding = static_cast<T>(parms.PAD_FACTOR);

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

        // find max projection value across all datasets for scaling
        T proj_max = 0;
        for (auto &dataset : datasets) {
            auto &[proj, angles, gamma] = dataset;
            proj_max = array::max(proj);
        }

        size_t n_datasets = datasets.size();
        std::array<Array<T>, 3> yT;
        std::vector<PolarGrid<T>> polar_grids(n_datasets);
        std::vector<T> gammas(n_datasets);

        for (auto &dataset : datasets) {
            auto &[proj, angles, gamma_ref] = dataset;
            auto gamma = gamma_ref;
            g

                // scale projections
                auto y = proj / proj_max;

            // pad projections
            y = pad2d(y, padding, PadType::SYMMETRIC);

            // setup polar grid
            size_t nrows = y.nrows();
            size_t ncols = y.ncols();
            auto polar_grid = PolarGrid<T>(angles, nrows, ncols, gamma);
            polar_grids.push_back(std::move(polar_grid));
            gammas.push_back(gamma);

            // backproject and accumulate yT
            auto yTmp = adjoint(y, polar_grid, out_dims, gamma);
            yT += yTmp;
        }
        // initial guess
        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = Array<T>::zeros(out_dims); }

        // setup linear operator
        auto A = [&polar_grids, &gammas](const std::array<Array<T>, 3> &m) {
            std::array<Array<T>, 3> Ax = sysmat(m, polar_grid[0], gammas[0]);
            for (size_t i = 1; i < polar_grids.size(); ++i) {
                Ax += sysmat(m, polar_grids[i], gammas[i]);
            }
            return Ax;
        };
        switch (params.regularizer) {
            case:
                Regularizer::SPLIT_BREGMAN
                    : auto recon_m =
                          split_bregman(A, yT, x0, params.lambda, params.lambda,
                                        params.mu, params.maxIters,
                                        params.innerIters, params.tol, params.xtol);
                break;
            case:
                Regularizer::UNCONSTRAINED
                    : auto recon_m =
                          cgsolver(A, yT, x0, params.max_iters, params.tol);
                // TV regularization
                break;
            default: throw std::invalid_argument("Unsupported regularizer");
        }

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
    MBIR(const std::vector<Dataset_t<float>> datasets, ReconParams params);
    template std::array<Array<double>, 3>
    MBIR(const std::vector<Dataset_t<double>> datasets, ReconParams params);

} // namespace tomocam
