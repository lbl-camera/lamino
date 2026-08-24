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
#include <cassert>
#include <format>
#include <functional>
#include <iostream>

#include "recon_params.h"

#include "array.h"
#include "array_ops.h"
#include "gpu/cufft_plan_cache.h"
#include "gpu/cufinufft_plan_cache.h"
#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/gpu_opt.h"
#include "gpu/padding.h"
#include "gpu/polar_grid.h"
#include "gpu/projection.h"
#include "gpu/vec_array.h"
#include "mask.h"

namespace tomocam::gpu {
    template <typename T>
    using Dataset_t = std::tuple<Array<T>, std::vector<T>, T>;

    template <typename T>
    std::array<Array<T>, 3> MBIR(const std::vector<Dataset_t<T>> &datasets,
                                 const ReconParams &params) {

        T proj_max = (T)0;
        for (auto &[projs, angles, gamma] : datasets) {
            proj_max = std::max(proj_max, tomocam::array::max(projs));
        }

        T padfac = static_cast<T>(params.PAD_FACTOR);

        // get recon dimensions from the params
        dims_t recon_dims = params.recon_dims;

        // extend recon dimensions to match padded projection size
        dims_t proj_dims = std::get<0>(datasets[0]).dims();
        dims_t out_dims = {recon_dims.n1, static_cast<size_t>(proj_dims.n2 * padfac),
                           static_cast<size_t>(proj_dims.n3 * padfac)};

        // move data back to host after all GPU work is done; declared outside the
        // device scope so it survives past cudaDeviceReset()
        std::array<Array<T>, 3> recon_host;
        {
            // setup the linear system for all datasets
            size_t n_datasets = datasets.size();
            std::vector<PolarGrid<T>> polar_grids(n_datasets);
            std::vector<T> gammas(n_datasets, (T)0);
            VecArray<T> yT{DeviceArray<T>(out_dims), DeviceArray<T>(out_dims),
                           DeviceArray<T>(out_dims)};

            for (size_t i = 0; i < n_datasets; ++i) {
                auto &[projs, angles, gamma_ref] = datasets[i];
                auto gamma = gamma_ref;
                gammas[i] = gamma;

                // move data to device and normalize
                DeviceArray<T> y(projs);
                y /= proj_max;

                // zero-pad projections by sqrt(2) to avoid aliasing
                y = pad2d(y, padfac, PadType::SYMMETRIC);

                // setup polar grid
                size_t nrows = y.nrows();
                size_t ncols = y.ncols();
                polar_grids[i] =
                    std::move(PolarGrid<T>(angles, gamma, nrows, ncols));

                // backproject y to get A^T y for optimization
                auto yTmp = adjoint(y, polar_grids[i], out_dims, gamma);
                for (size_t j = 0; j < 3; ++j) { yT[j] += yTmp[j]; }
            }

            // setup the linear operator for all datasets
            opt::gpuFunction<T> A = [&polar_grids,
                                     &gammas](const gpu::VecArray<T> &x) {
                auto Ax = sysmat<T>(x, polar_grids[0], gammas[0]);
                for (size_t i = 1; i < gammas.size(); ++i) {
                    auto tmp = sysmat<T>(x, polar_grids[i], gammas[i]);
                    for (size_t j = 0; j < 3; ++j) { Ax[j] += tmp[j]; }
                }
                return Ax;
            };

            // initialize solution with backprojection of yT
            VecArray<T> x0{DeviceArray<T>(out_dims), DeviceArray<T>(out_dims),
                           DeviceArray<T>(out_dims)};

            VecArray<T> recon;
            switch (params.regularizer) {
                case Regularizer::UNCONSTRAINED: {
                    std::cout << "Starting unconstrained reconstruction with CG on "
                                 "GPU ...\n";
                    recon = opt::cgsolver<T>(A, yT, x0, params.maxIters, params.tol,
                                             params.xtol, out_dims);
                    break;
                }
                case Regularizer::SPLIT_BREGMAN: {
                    std::cout
                        << "Starting MBIR with Split-Bregman method on GPU ...\n";
                    recon = opt::split_bregman<T>(
                        A, yT, x0, params.lambda, params.mu, params.maxIters,
                        params.innerIters, params.tol, params.xtol, out_dims);

                    break;
                }
                default: throw std::invalid_argument("Unsupported optimizer type");
            }
            // crop to original dimensions
            for (size_t i = 0; i < 3; ++i) {
                recon[i] = crop3d(recon[i], recon_dims, PadType::SYMMETRIC);
            }

            for (size_t i = 0; i < 3; ++i) { recon_host[i] = recon[i].to_host(); }

            // cleanup plan caches; all device objects destroyed at end of this block
            // before cudaDeviceReset() is called below
            tomocam::gpu::nufft::plans::cache<float>.clear();
            tomocam::gpu::fft::plans::cache<float>.clear();
        }
        cudaDeviceReset();

        return recon_host;
    }

    // explicit template instantiation
    template std::array<Array<float>, 3>
    MBIR(const std::vector<Dataset_t<float>> &datasets, const ReconParams &params);
    template std::array<Array<double>, 3>
    MBIR(const std::vector<Dataset_t<double>> &datasets, const ReconParams &params);

} // namespace tomocam::gpu
