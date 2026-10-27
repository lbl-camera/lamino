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
#include "logger.h"
#include "mask.h"
#include "optimize.h"
#include "padding.h"
#include "polar_grid.h"
#include "projection.h"
#include "recon_params.h"
#include "timer.h"
#include "toeplitz.h"

namespace tomocam {

    template <typename T>
    std::array<Array<T>, 3> MBIR(const std::vector<Dataset_t<T>> &datasets,
                                 const ReconParams &params, const Array<T> &mask) {

        // create logger (defaults to stdout)
        Logger logger(params.logMode, params.logfile);

        // padding factor
        T padfac = static_cast<T>(params.PAD_FACTOR);

        dims_t proj_dims = datasets[0].projs.dims();
        dims_t output_dims = params.recon_dims;
        // The detector/PolarGrid is padded to P = max(recon n2/n3, sqrt(2) *
        // detector); the volume is V = max(odd detector, recon n2/n3) <= P
        // (mirrors src/gpu/mbir.cu).
        size_t grid_n2 = recon_grid_dim<T>(proj_dims.n2, output_dims.n2, padfac);
        size_t grid_n3 = recon_grid_dim<T>(proj_dims.n3, output_dims.n3, padfac);
        dims_t recon_dims = {padded_dim(output_dims.n1, padfac),
                             recon_vol_dim(proj_dims.n2, output_dims.n2),
                             recon_vol_dim(proj_dims.n3, output_dims.n3)};

        std::array<Array<T>, 3> yT;
        for (size_t i = 0; i < 3; ++i) { yT[i] = Array<T>::zeros(recon_dims); }

        size_t n_datasets = datasets.size();
        std::vector<PolarGrid<T>> polar_grids(n_datasets);

        // start preprocessing timer
        Timer t0;
        t0.start();
        T proj_max = 0.0;
        for (const auto &ds : datasets) {
            proj_max = std::max(proj_max, array::max(ds.projs));
        }

        for (size_t j = 0; j < n_datasets; ++j) {
            const auto &ds = datasets[j];

            auto y = ds.projs / proj_max;
            y = pad2d(y, grid_n2, grid_n3, PadType::SYMMETRIC);

            size_t nrows = y.nrows();
            size_t ncols = y.ncols();
            polar_grids[j] =
                std::move(PolarGrid<T>(ds.angles, nrows, ncols, ds.gamma, ds.beta));

            auto yTmp =
                adjoint(y, polar_grids[j], recon_dims, ds.gamma, ds.beta, ds.shifts);
            for (size_t i = 0; i < 3; ++i) { yT[i] += yTmp[i]; }
        }

        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = yT[i].clone(); }

        // Precompute the Toeplitz PSF kernels once per dataset grid, so A^T A
        // becomes a set of FFT convolutions instead of a NUFFT
        // forward+adjoint pair per solver iteration.
        std::vector<cpu::ToeplitzVectorOp<T>> toeplitz_ops;
        toeplitz_ops.reserve(n_datasets);
        for (size_t j = 0; j < n_datasets; ++j) {
            toeplitz_ops.emplace_back(polar_grids[j], recon_dims);
        }

        opt::Function<T> A = [&toeplitz_ops](const std::array<Array<T>, 3> &m) {
            std::array<Array<T>, 3> Ax = cpu::sysmat(m, toeplitz_ops[0]);
            for (size_t i = 1; i < toeplitz_ops.size(); ++i) {
                auto tmp = cpu::sysmat(m, toeplitz_ops[i]);
                for (size_t j = 0; j < 3; ++j) { Ax[j] += tmp[j]; }
            }
            return Ax;
        };

        t0.stop();
        logger.log(std::format("Preprocessing took {} seconds", t0.seconds()));

        // the const-mag support, centred in the (padded) solver volume; the
        // padding stays outside the support and is left free
        Array<T> solver_mask = embed_mask(mask, recon_dims);
        // the mask only selects where |m| = M0 is enforced: voxels outside it
        // are free, and are forced to 0 only if zero_outside is set (a mask
        // file alone does not zero them)
        const bool use_support = params.zero_outside;
        const Array<T> no_support;
        const Array<T> &support = use_support ? solver_mask : no_support;

        std::array<Array<T>, 3> recon_m;
        switch (params.regularizer) {
            case Regularizer::SPLIT_BREGMAN:
                recon_m = opt::split_bregman<T>(
                    A, yT, x0, params.lambda, params.mu, params.maxIters,
                    params.innerIters, params.tol, params.xtol, support, &logger);
                break;
            case Regularizer::CONST_MAG:
                recon_m = opt::constmag_l2<T>(A, yT, x0, params.lambda_l2, params.mu,
                                              params.maxIters, params.innerIters,
                                              params.tol, params.xtol, solver_mask,
                                              &logger, use_support);
                break;
            case Regularizer::COMBINED: {
                opt::RegWeights<T> w{params.const_mag, params.lambda,
                                     params.lambda_l2};
                opt::MagOptions<T> mag;
                mag.rho = params.rho;
                mag.rho_growth = params.rho_growth;
                mag.rho_max = params.rho_max;
                mag.Ms = params.Ms;
                mag.estimate_Ms = params.estimate_Ms;
                mag.zero_outside = use_support;
                recon_m = opt::regularized_admm<T>(
                    A, yT, x0, w, params.mu, params.maxIters, params.innerIters,
                    params.tol, params.xtol, solver_mask, &logger, mag);
                break;
            }
            case Regularizer::UNCONSTRAINED:
                recon_m = opt::cgsolver<T>(A, yT, x0, params.maxIters, params.tol,
                                           params.xtol, support, T(0), &logger);
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
    MBIR(const std::vector<Dataset_t<float>> &datasets, const ReconParams &params,
         const Array<float> &mask);
    template std::array<Array<double>, 3>
    MBIR(const std::vector<Dataset_t<double>> &datasets, const ReconParams &params,
         const Array<double> &mask);

    /*--------------------------------------------------------------------------------
     * MBIR2
     * --------------------------------------------------------------------------------*/
    template <typename T>
    std::array<Array<T>, 3> MBIR2(const std::vector<Dataset_t<T>> &datasets,
                                  const ReconParams &params, const Array<T> &mask) {

        Logger logger(params.logMode, params.logfile);

        T padfac = static_cast<T>(params.PAD_FACTOR);

        dims_t proj_dims = datasets[0].projs.dims();
        dims_t output_dims = params.recon_dims;
        // The detector/PolarGrid is padded to P = max(recon n2/n3, sqrt(2) *
        // detector); the volume is V = max(odd detector, recon n2/n3) <= P
        // (mirrors src/gpu/mbir.cu).
        size_t grid_n2 = recon_grid_dim<T>(proj_dims.n2, output_dims.n2, padfac);
        size_t grid_n3 = recon_grid_dim<T>(proj_dims.n3, output_dims.n3, padfac);
        dims_t recon_dims = {padded_dim(output_dims.n1, padfac),
                             recon_vol_dim(proj_dims.n2, output_dims.n2),
                             recon_vol_dim(proj_dims.n3, output_dims.n3)};

        for (const auto &ds : datasets) {
            if (ds.projs.nrows() != proj_dims.n2 ||
                ds.projs.ncols() != proj_dims.n3) {
                throw std::invalid_argument(
                    "MBIR: all datasets must have the same projection size");
            }
        }
        logger.log(std::format("volume {}x{}x{}, projection/polar grid {}x{}\n",
                               recon_dims.n1, recon_dims.n2, recon_dims.n3, grid_n2,
                               grid_n3));

        T proj_max = 0.0;
        for (const auto &ds : datasets) {
            proj_max = std::max(proj_max, array::max(ds.projs));
        }

        std::vector<std::tuple<std::vector<T>, T, T>> angle_gamma_beta;
        size_t total_nangles = 0;
        size_t nrows = 0, ncols = 0;

        std::vector<Array<T>> padded;
        std::vector<std::array<T, 2>> shifts; // stacked in the same order as y
        bool has_shifts = false;
        for (const auto &ds : datasets) {
            auto y =
                pad2d(ds.projs / proj_max, grid_n2, grid_n3, PadType::SYMMETRIC);
            assert(nrows == 0 || (y.nrows() == nrows && y.ncols() == ncols));
            nrows = y.nrows();
            ncols = y.ncols();
            angle_gamma_beta.push_back({ds.angles, ds.gamma, ds.beta});
            total_nangles += ds.angles.size();
            padded.push_back(std::move(y));

            if (ds.shifts.empty()) {
                shifts.insert(shifts.end(), ds.angles.size(), {T(0), T(0)});
            } else {
                assert(ds.shifts.size() == ds.angles.size());
                shifts.insert(shifts.end(), ds.shifts.begin(), ds.shifts.end());
                for (const auto &s : ds.shifts) {
                    if (s[0] != T(0) || s[1] != T(0)) has_shifts = true;
                }
            }
        }

        Array<T> y_stacked(dims_t{total_nangles, nrows, ncols});
        size_t offset = 0;
        for (size_t j = 0; j < padded.size(); ++j) {
            size_t n = padded[j].nslices();
            auto src = padded[j].slice(0, n);
            auto dst = y_stacked.slice(offset, offset + n);
            std::copy(std::execution::par_unseq, src.begin(), src.end(),
                      dst.begin());
            offset += n;
        }

        PolarGrid<T> pg(angle_gamma_beta, nrows, ncols);

        // Backproject once with the unified grid
        // skip the phase ramp entirely when no dataset has a nonzero shift
        if (!has_shifts) shifts.clear();
        auto yT = adjoint(y_stacked, pg, recon_dims, shifts);

        // initial guess
        std::array<Array<T>, 3> x0;
        for (size_t i = 0; i < 3; ++i) { x0[i] = yT[i].clone(); }

        // Precompute the Toeplitz PSF kernels once for the unified grid
        cpu::ToeplitzVectorOp<T> toeplitz_op(pg, recon_dims);

        opt::Function<T> A = [&toeplitz_op](const std::array<Array<T>, 3> &m) {
            return cpu::sysmat(m, toeplitz_op);
        };

        // the const-mag support, centred in the (padded) solver volume; the
        // padding stays outside the support and is left free
        Array<T> solver_mask = embed_mask(mask, recon_dims);
        // the mask only selects where |m| = M0 is enforced: voxels outside it
        // are free, and are forced to 0 only if zero_outside is set (a mask
        // file alone does not zero them)
        const bool use_support = params.zero_outside;
        const Array<T> no_support;
        const Array<T> &support = use_support ? solver_mask : no_support;

        std::array<Array<T>, 3> recon_m;
        switch (params.regularizer) {
            case Regularizer::SPLIT_BREGMAN:
                recon_m = opt::split_bregman<T>(
                    A, yT, x0, params.lambda, params.mu, params.maxIters,
                    params.innerIters, params.tol, params.xtol, support, &logger);
                break;
            case Regularizer::CONST_MAG:
                // constrain only the block that survives the crop below; the
                // padding stays free to absorb out-of-view contributions
                recon_m = opt::constmag_l2<T>(A, yT, x0, params.lambda_l2, params.mu,
                                              params.maxIters, params.innerIters,
                                              params.tol, params.xtol, solver_mask,
                                              &logger, use_support);
                break;
            case Regularizer::COMBINED: {
                opt::RegWeights<T> w{params.const_mag, params.lambda,
                                     params.lambda_l2};
                opt::MagOptions<T> mag;
                mag.rho = params.rho;
                mag.rho_growth = params.rho_growth;
                mag.rho_max = params.rho_max;
                mag.Ms = params.Ms;
                mag.estimate_Ms = params.estimate_Ms;
                mag.zero_outside = use_support;
                recon_m = opt::regularized_admm<T>(
                    A, yT, x0, w, params.mu, params.maxIters, params.innerIters,
                    params.tol, params.xtol, solver_mask, &logger, mag);
                break;
            }
            case Regularizer::UNCONSTRAINED:
                recon_m = opt::cgsolver<T>(A, yT, x0, params.maxIters, params.tol,
                                           params.xtol, support, T(0), &logger);
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
    MBIR2(const std::vector<Dataset_t<float>> &datasets, const ReconParams &params,
          const Array<float> &mask);
    template std::array<Array<double>, 3>
    MBIR2(const std::vector<Dataset_t<double>> &datasets, const ReconParams &params,
          const Array<double> &mask);

} // namespace tomocam
