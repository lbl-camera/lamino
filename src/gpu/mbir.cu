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
#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/gpu_opt.h"
#include "gpu/padding.h"
#include "gpu/polar_grid.h"
#include "gpu/projection.h"

namespace tomocam::gpu {
    template <typename T>
    using Dataset_t = std::tuple<Array<T>, std::vector<T>, T>;

    template <typename T>
    std::array<Array<T>, 3> MBIR(const std::vector<Dataset_t<T>> &datasets,
                                       const ReconParams &params) {

        T scale = (T)0;
        for (auto &[projs, angles, gamma, recon_dims] : datasets) {
            scale = std::max(scale, array::max(projs));
        }

        // extend recon dimenstion by PAD_FACTOR
        T padding = static_cast<T>(params.PAD_FACTOR) - (T)1.0;
        dims_t out_dims = params.recon_dims;
        size_t n1_pad = 2 * (static_cast<size_t>(recon_dims.n1 * padding) / 2);
        out_dims.n1 += n1_pad;
        size_t n2_pad = 2 * (static_cast<size_t>(recon_dims.n2 * padding) / 2);
        out_dims.n2 += n2_pad;
        size_t n3_pad = 2 * (static_cast<size_t>(recon_dims.n3 * padding) / 2);
        out_dims.n3 += n3_pad;

        // setup the linear system for all datasets
        size_t n_datasets = datasets.size();
        std::vector<PolarGrid<T>> polar_grids;
        std::vector<T> gammas(n_datasets);
        std::array<DeviceArray<T>, 3> yT;

        for (size_t i = 0; i < n_datasets; ++i) {
            auto &[projs, angles, gamma_ref] = datasets[i];
            T gammas[i] = gramma_ref;

            // normalize projections
            DeviceArray<T> y(projs);
            y =/ scale;

            // zero-pad projections by sqrt(2) to avoid aliasing
            float padding = static_cast<T>(params.PAD_FACTOR);
            y = pad2d(y, padding, PadType::SYMMETRIC);

            // setup polar grid
            size_t nrows = y.nrows();
            size_t ncols = y.ncols();
            polar_grids[i] =
                std::move(PolarGrid<T>(angles, gammas[i], nrows, ncols));

            // backproject y to get A^T y for optimization
            auto yTmp = backproj(y, polar_grids[i], out_dims);
            for (size_t j = 0; j < 3; ++j) { yT[j] += yTmp[j]; }
        }

        // setup the linear operator for all datasets
        opt::gpuFunction<T> A = [&](const std::array<DeviceArray<T>, 3> &x) {
            auto Ax = sysmat<T>(x, polar_grids[0], gammas[0]);
            for (size_t i = 1; i < n_datasets; ++i) {
                Ax += sysmat<T>(x, polar_grids[i], gammas[i]);
            }
            return Ax;
        };

        // initialize solution with zeros
        auto x0 = std::array<DeviceArray<T>, 3>{DeviceArray<T>(out_dims),
                                                DeviceArray<T>(out_dims),
                                                DeviceArray<T>(out_dims)};

        // run optimization
        auto recon = std::array<DeviceArray<T>, 3>{DeviceArray<T>(out_dims),
                                                   DeviceArray<T>(out_dims),
                                                   DeviceArray<T>(out_dims)};
        switch (params.regularizer) {
            case Regularizer::UNCONSTRAINED: {
                std::cout
                    << "Starting unconstrained reconstruction with CG on GPU ...\n";
                recon = opt::cgsolver<T>(A, yT, x0, params.maxIters, params.tol,
                                         params.xtol);
                break;
            }
            case Regularizer::SPLIT_BREGMAN: {
                std::cout << "Starting MBIR with Split-Bregman method on GPU ...\n";
                recon = opt::split_bregman<T>(A, yT, x0, params.lambda, params.mu,
                                              params.maxIters, params.innerIters,
                                              params.tol, params.xtol);

                break;
            }
            default: throw std::invalid_argument("Unsupported optimizer type");
        }

        // crop to original dimensions
        for (size_t i = 0; i < 3; ++i) {
            recon[i] = crop3d(recon[i], recon_dims, PadType::SYMMETRIC);
        }
    }

    // explicit template instantiation
    template DeviceArray<float> MBIR<float>(const DeviceArray<float> &projs,
                                            const std::vector<float> &angles,
                                            float gamma, const dims_t &recon_dims,
                                            const ReconParams &cfg);
    template DeviceArray<double> MBIR<double>(const DeviceArray<double> &projs,
                                              const std::vector<double> &angles,
                                              double gamma, const dims_t &recon_dims,
                                              const ReconParams &cfg);

} // namespace tomocam::gpu
