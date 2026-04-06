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
#include "mask.h"
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
    std::array<Array<T>, 3> MBIR1(const std::vector<Dataset_t<T>> &datasets,
                                  const dims_t &recon_dims,
                                  const ReconParams &recon_params) {

        // padding factor
        T padding = static_cast<T>(recon_params.PAD_FACTOR - 1.0);

        // adjust reconstruction dimensions
        dims_t out_dims = recon_dims;
        size_t n1_pad = 2 * (static_cast<size_t>(recon_dims.n1 * padding) / 2);
        out_dims.n1 += n1_pad;

        size_t n2_pad = 2 * (static_cast<size_t>(recon_dims.n2 * padding) / 2);
        out_dims.n2 += n2_pad;

        size_t n3_pad = 2 * (static_cast<size_t>(recon_dims.n3 * padding) / 2);
        out_dims.n3 += n3_pad;

        // setup system matrices and backprojections
        size_t n_datasets = datasets.size();
        std::vector<PolarGrid<T>> polar_grids(n_datasets);

        // array for backprojected measurements
        std::array<Array<T>, 3> yT;
        for (size_t i = 0; i < 3; ++i) { yT[i] = Array<T>::zeros(out_dims); }
        std::vector<T> gammas(n_datasets);
        for (size_t j = 0; j < n_datasets; ++j) {
            auto &[proj, angles, gamma_ref] = datasets[j];
            gammas[j] = gamma_ref;

            // remove nans and infs
            auto proj2 = mask_infs_nans(proj);

            // normalize projections
            T proj_max = array::max(proj2);
            auto y = proj2 / proj_max;

            // zero-pad projections by sqrt(2) to avoid aliasing
            y = pad2d(y, padding, PadType::SYMMETRIC);

            // setup polar grid
            size_t nrows = y.nrows();
            size_t ncols = y.ncols();
            polar_grids[j] =
                std::move(PolarGrid<T>(angles, nrows, ncols, gammas[j]));

            // backproject measurements to get yT
            auto yTmp = adjoint(y, polar_grids[j], out_dims, gammas[j]);
            for (size_t i = 0; i < 3; ++i) { yT[i] += yTmp[i]; }
        }

        // setup the linear system for CG solver
        opt::Function<T> A = [&polar_grids,
                              &gammas](const std::array<Array<T>, 3> &x) {
            std::array<Array<T>, 3> Ax = sysmat(x, polar_grids[0], gammas[0]);
            for (size_t j = 1; j < polar_grids.size(); ++j) {
                auto Atmp = sysmat(x, polar_grids[j], gammas[j]);
                for (size_t i = 0; i < 3; ++i) { Ax[i] += Atmp[i]; }
            }
            return Ax;
        };

        // initial guess
        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = Array<T>::zeros(out_dims); }

        // demagnetization constraint weight
        // Lambda controls divergence-free constraint: higher values enforce ∇·M ≈ 0
        // Typical range: 0.001 - 0.1 (relative to data fidelity term)
        T lambda = recon_params.lambda;

        // solve linear system using CG solver with demagnetization constraint
        std::array<Array<T>, 3> recon_m = opt::cgsolver<T>(
            A, yT, x0, recon_params.maxIters, recon_params.tol, lambda);

        // crop to original dimensions
        std::array<Array<T>, 3> recon_magnetisation;
        for (size_t i = 0; i < 3; ++i) {
            recon_magnetisation[i] =
                crop3d(recon_m[i], recon_dims, PadType::SYMMETRIC);
        }

        // transpose to match expected output format
        for (size_t i = 0; i < 3; ++i) {
            recon_magnetisation[i] =
                array::transpose(recon_magnetisation[i], {1, 2, 0});
        }
        return recon_magnetisation;
    }

    // Explicit template instantiations
    template std::array<Array<float>, 3>
    MBIR1(const std::vector<Dataset_t<float>> &datasets, const dims_t &recon_dims,
          const ReconParams &recon_params);
    template std::array<Array<double>, 3>
    MBIR1(const std::vector<Dataset_t<double>> &datasets, const dims_t &recon_dims,
          const ReconParams &recon_params);
} // namespace tomocam
