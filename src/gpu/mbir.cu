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
#include <tuple>
#include <vector>

#include "recon_params.h"

#include "array.h"
#include "array_ops.h"
#include "logger.h"
#include "mask.h"

#include "gpu/cufft_plan_cache.h"
#include "gpu/cufinufft_plan_cache.h"
#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/gpu_opt.h"
#include "gpu/padding.h"
#include "gpu/polar_grid.h"
#include "gpu/projection.h"
#include "gpu/toeplitz.h"
#include "gpu/tomocam.h"
#include "gpu/vec_array.h"

namespace tomocam::gpu {

    template <typename T>
    std::array<Array<T>, 3> MBIR(const std::vector<Dataset_t<T>> &datasets,
                                 const ReconParams &params) {

        // create logger
        Logger logger(params.logMode, params.logfile);

        T proj_max = (T)0;
        for (const auto &ds : datasets) {
            proj_max = std::max(proj_max, tomocam::array::max(ds.projs));
        }

        // pad solution and projection data
        T padfac = static_cast<T>(params.PAD_FACTOR);
        dims_t recon_dims = params.recon_dims;

        // extend recon dimensions: n1 padded by PAD_FACTOR (prevents NUFFT
        // z-aliasing), n2/n3 derived from padded projection dimensions
        dims_t proj_dims = datasets[0].projs.dims();
        // The detector/PolarGrid stays padded (avoids 2D FFT aliasing), but the
        // volume in n2/n3 stays at the unpadded projection size: the type-1
        // NUFFT evaluates the backprojection and the Toeplitz PSF exactly on
        // any mode window, so the PSF grid is 2*n-1 of the unpadded size
        // instead of 2*(padfac*n)-1 (~2x fewer voxels, less cufinufft memory).
        dims_t out_dims = {padded_dim<T>(recon_dims.n1, padfac), proj_dims.n2,
                           proj_dims.n3};

        // move data back to host after all GPU work is done; declared outside the
        // device scope so it survives past cudaDeviceReset()
        std::array<Array<T>, 3> recon_host;
        {
            // stack all datasets into one unified polar grid (mirrors CPU
            // MBIR2): a single backprojection and a single Toeplitz PSF
            // instead of one per dataset.
            size_t n_datasets = datasets.size();
            std::vector<std::tuple<std::vector<T>, T, T>> angle_gamma_beta;
            size_t total_nangles = 0;
            for (const auto &ds : datasets) {
                angle_gamma_beta.push_back({ds.angles, ds.gamma, ds.beta});
                total_nangles += ds.angles.size();
            }

            // per-projection center-of-rotation shifts, stacked in the same
            // order as y; skipped entirely when no dataset has a nonzero shift
            std::vector<T> hx(total_nangles, T(0)), hy(total_nangles, T(0));
            bool has_shifts = false;

            DeviceArray<T> y_stacked;
            size_t nrows = 0, ncols = 0, offset = 0;
            for (size_t i = 0; i < n_datasets; ++i) {
                const auto &ds = datasets[i];

                // move data to device and normalize
                DeviceArray<T> y(ds.projs);
                y /= proj_max;

                // zero-pad projections by sqrt(2) to avoid aliasing
                y = pad2d(y, padfac, PadType::SYMMETRIC);

                if (i == 0) {
                    nrows = y.nrows();
                    ncols = y.ncols();
                    y_stacked = DeviceArray<T>(dims_t{total_nangles, nrows, ncols});
                }
                assert(y.nrows() == nrows && y.ncols() == ncols);
                copyD2D(y_stacked.data() + offset * nrows * ncols, y.data(),
                        y.bytes());

                if (!ds.shifts.empty()) {
                    assert(ds.shifts.size() == ds.angles.size());
                    for (size_t k = 0; k < ds.shifts.size(); ++k) {
                        hx[offset + k] = ds.shifts[k][0];
                        hy[offset + k] = ds.shifts[k][1];
                        if (hx[offset + k] != T(0) || hy[offset + k] != T(0)) {
                            has_shifts = true;
                        }
                    }
                }
                offset += y.nslices();
            }

            gpu::PolarGrid<T> pg(angle_gamma_beta, nrows, ncols);

            thrust::device_vector<T> shift_dx, shift_dy;
            if (has_shifts) {
                shift_dx = thrust::device_vector<T>(hx);
                shift_dy = thrust::device_vector<T>(hy);
            }

            // backproject once with the unified grid to get A^T y
            VecArray<T> yT = adjoint(y_stacked, pg, out_dims, shift_dx, shift_dy);
            y_stacked = DeviceArray<T>();

            // Precompute the Toeplitz PSF kernels once for the unified grid, so
            // A^T A becomes a set of FFT convolutions instead of a NUFFT
            // forward+adjoint pair per solver iteration.
            gpu::ToeplitzVectorOp<T> toeplitz_op(pg, out_dims,
                                                 gpu::ToeplitzMode::SEQUENTIAL);

            // setup the linear operator for all datasets
            opt::gpuFunction<T> A = [&toeplitz_op](const gpu::VecArray<T> &x) {
                return sysmat<T>(x, toeplitz_op);
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
                                             params.xtol, out_dims, T(0), &logger);
                    break;
                }
                case Regularizer::SPLIT_BREGMAN: {
                    std::cout
                        << "Starting MBIR with Split-Bregman method on GPU ...\n";
                    recon = opt::split_bregman<T>(A, yT, x0, params.lambda,
                                                  params.mu, params.maxIters,
                                                  params.innerIters, params.tol,
                                                  params.xtol, out_dims, &logger);

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
