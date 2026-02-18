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
    std::array<Array<T>, 3> MBIR(const std::vector<Dataset_t<T>> &datasets,
                                 const dims_t &recon_dims,
                                 const ReconParams &recon_params) {

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

        // setup system matrices and backprojections
        size_t n_datasets = datasets.size();
        std::vector<opt::Function<T>> sysmats;
        std::vector<std::array<Array<T>, 3>> yTs;
        std::vector<T> gammas;
        for (auto &dataset : datasets) {
            auto &[proj, angles, gamma_ref] = dataset;
            T gamma = gamma_ref;
            gammas.push_back(gamma);

            // normalize projections
            T proj_max = array::max(proj);
            auto y = proj / proj_max;

            // zero-pad projections by sqrt(2) to avoid aliasing
            y = pad2d(y, padding, PadType::SYMMETRIC);

            // setup polar grid
            size_t nrows = y.nrows();
            size_t ncols = y.ncols();
            auto polar_grid =
                std::make_shared<PolarGrid<T>>(angles, nrows, ncols, gamma);

            // backproject measurements to get yT
            auto yT = adjoint(y, *polar_grid.get(), out_dims, gamma);
            yTs.push_back(std::move(yT));

            // setup gradient operator
            opt::Function<T> A = [pg = polar_grid,
                                  gamma](const std::array<Array<T>, 3> &x) {
                // gradient data
                return sysmat(x, *pg.get(), gamma);
            };
            sysmats.push_back(std::move(A));
        }

        // initial guess
        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = Array<T>::ones(out_dims) * 0.9; }
        auto recon_m = opt::split_bregman<T>(sysmats, yTs, x0, recon_params.lambda,
                                             recon_params.mu, recon_params.maxIters,
                                             recon_params.innerIters,
                                             recon_params.tol, recon_params.xtol);

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
    MBIR(const std::vector<Dataset_t<float>> &datasets, const dims_t &recon_dims,
         const ReconParams &recon_params);
    template std::array<Array<double>, 3>
    MBIR(const std::vector<Dataset_t<double>> &datasets, const dims_t &recon_dims,
         const ReconParams &recon_params);

} // namespace tomocam
