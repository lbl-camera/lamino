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
#include <execution>
#include <format>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "mask.h"
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

        dims_t proj_dims = datasets[0].projs.dims();
        dims_t output_dims = params.recon_dims;
        dims_t recon_dims = {output_dims.n1,
                             static_cast<size_t>(proj_dims.n2 * padfac),
                             static_cast<size_t>(proj_dims.n3 * padfac)};

        std::array<Array<T>, 3> yT;
        for (size_t i = 0; i < 3; ++i) { yT[i] = Array<T>::zeros(recon_dims); }

        size_t n_datasets = datasets.size();
        std::vector<PolarGrid<T>> polar_grids(n_datasets);
        std::vector<T> betas(n_datasets);

        T proj_max = 0.0;
        for (const auto &ds : datasets) {
            proj_max = std::max(proj_max, array::max(ds.projs));
        }

        for (size_t j = 0; j < n_datasets; ++j) {
            const auto &ds = datasets[j];
            betas[j] = ds.beta;

            auto y = ds.projs / proj_max;
            y = pad2d(y, padfac, PadType::SYMMETRIC);

            size_t nrows = y.nrows();
            size_t ncols = y.ncols();
            polar_grids[j] =
                std::move(PolarGrid<T>(ds.angles, nrows, ncols, ds.gamma, ds.beta));

            auto yTmp = adjoint(y, polar_grids[j], recon_dims, ds.beta, ds.shifts);
            for (size_t i = 0; i < 3; ++i) { yT[i] += yTmp[i]; }
        }

        {
            auto supp_mask = mask_support<T>(recon_dims, output_dims);
            for (size_t i = 0; i < 3; ++i) { yT[i] *= supp_mask; }
        }

        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = yT[i].clone(); }

        opt::Function<T> A = [&polar_grids,
                              &betas](const std::array<Array<T>, 3> &m) {
            std::array<Array<T>, 3> Ax = sysmat(m, polar_grids[0], betas[0]);
            for (size_t i = 1; i < polar_grids.size(); ++i) {
                auto tmp = sysmat(m, polar_grids[i], betas[i]);
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

    template <typename T>
    std::array<Array<T>, 3> MBIR2(const std::vector<Dataset_t<T>> &datasets,
                                   const ReconParams &params) {

        T padfac = static_cast<T>(params.PAD_FACTOR);

        dims_t proj_dims = datasets[0].projs.dims();
        dims_t output_dims = params.recon_dims;
        dims_t recon_dims = {output_dims.n1,
                             static_cast<size_t>(proj_dims.n2 * padfac),
                             static_cast<size_t>(proj_dims.n3 * padfac)};

        T proj_max = 0.0;
        for (const auto &ds : datasets) {
            proj_max = std::max(proj_max, array::max(ds.projs));
        }

        std::vector<std::tuple<std::vector<T>, T, T>> angle_gamma_beta;
        size_t total_nangles = 0;
        size_t nrows = 0, ncols = 0;

        std::vector<Array<T>> padded;
        for (const auto &ds : datasets) {
            auto y = pad2d(ds.projs / proj_max, padfac, PadType::SYMMETRIC);
            assert(nrows == 0 || (y.nrows() == nrows && y.ncols() == ncols));
            nrows = y.nrows();
            ncols = y.ncols();
            angle_gamma_beta.push_back({ds.angles, ds.gamma, ds.beta});
            total_nangles += ds.angles.size();
            padded.push_back(std::move(y));
        }

        Array<T> y_stacked(dims_t{total_nangles, nrows, ncols});
        size_t offset = 0;
        for (size_t j = 0; j < padded.size(); ++j) {
            size_t n = padded[j].nslices();
            auto src = padded[j].slice(0, n);
            auto dst = y_stacked.slice(offset, offset + n);
            std::copy(std::execution::par_unseq, src.begin(), src.end(), dst.begin());
            offset += n;
        }

        PolarGrid<T> pg(angle_gamma_beta, nrows, ncols);

        // Backproject once with the unified grid
        auto yT = adjoint(y_stacked, pg, recon_dims);

        // Zero yT outside the original (non-padded) support region
        {
            auto supp_mask = mask_support<T>(recon_dims, output_dims);
            for (size_t i = 0; i < 3; ++i) { yT[i] *= supp_mask; }
        }

        // Initial guess
        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = yT[i].clone(); }

        // System matrix: one call to the per-angle-gamma overload
        opt::Function<T> A = [&pg](const std::array<Array<T>, 3> &m) {
            return sysmat(m, pg);
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
                break;
            default: throw std::invalid_argument("Unsupported regularizer");
        }

        std::array<Array<T>, 3> recon_magnetisation;
        for (size_t i = 0; i < 3; ++i) {
            recon_magnetisation[i] =
                crop3d(recon_m[i], output_dims, PadType::SYMMETRIC);
        }
        return recon_magnetisation;
    }

    template std::array<Array<float>, 3>
    MBIR2(const std::vector<Dataset_t<float>> &datasets, const ReconParams &params);
    template std::array<Array<double>, 3>
    MBIR2(const std::vector<Dataset_t<double>> &datasets, const ReconParams &params);

} // namespace tomocam
