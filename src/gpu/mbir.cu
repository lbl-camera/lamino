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
                                 const ReconParams &params, const Array<T> &mask) {

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
        // The detector/PolarGrid is padded (avoids 2D FFT aliasing) to
        // P = max(recon n2/n3, sqrt(2) * detector), so recon_dims may exceed the
        // detector. The volume in n2/n3 is V = max(odd detector, recon n2/n3),
        // V <= P: the type-1 NUFFT evaluates the backprojection and the
        // Toeplitz PSF exactly on any mode window, so the PSF grid is 2*V-1
        // instead of 2*P-1 (less cufinufft memory).
        for (const auto &ds : datasets) {
            if (ds.projs.nrows() != proj_dims.n2 ||
                ds.projs.ncols() != proj_dims.n3) {
                throw std::invalid_argument(
                    "MBIR: all datasets must have the same projection size");
            }
        }
        size_t grid_n2 = recon_grid_dim<T>(proj_dims.n2, recon_dims.n2, padfac);
        size_t grid_n3 = recon_grid_dim<T>(proj_dims.n3, recon_dims.n3, padfac);
        dims_t out_dims = {padded_dim<T>(recon_dims.n1, padfac),
                           recon_vol_dim(proj_dims.n2, recon_dims.n2),
                           recon_vol_dim(proj_dims.n3, recon_dims.n3)};
        logger.log(std::format("volume {}x{}x{}, projection/polar grid {}x{}\n",
                               out_dims.n1, out_dims.n2, out_dims.n3, grid_n2,
                               grid_n3));

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

                // zero-pad projections to the polar grid size to avoid aliasing
                y = pad2d(y, grid_n2, grid_n3, PadType::SYMMETRIC);

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
            VecArray<T> x0 = yT.clone();

            // const-mag support, centred in the padded volume; the padding
            // stays outside the support and is left free
            DeviceArray<T> solver_mask(embed_mask(mask, out_dims));
            // the mask only selects where |m| = M0 is enforced: voxels outside it
            // are free, and are forced to 0 only if zero_outside is set (a mask
            // file alone does not zero them)
            const bool use_support = params.zero_outside;
            const DeviceArray<T> no_support;
            const DeviceArray<T> &support = use_support ? solver_mask : no_support;

            VecArray<T> recon;
            switch (params.regularizer) {
                case Regularizer::UNCONSTRAINED: {
                    std::cout << "Starting unconstrained reconstruction with CG on "
                                 "GPU ...\n";
                    recon = opt::cgsolver<T>(A, yT, x0, params.maxIters, params.tol,
                                             params.xtol, support, T(0), &logger);
                    break;
                }
                case Regularizer::SPLIT_BREGMAN: {
                    std::cout
                        << "Starting MBIR with Split-Bregman method on GPU ...\n";
                    recon = opt::split_bregman<T>(A, yT, x0, params.lambda,
                                                  params.mu, params.maxIters,
                                                  params.innerIters, params.tol,
                                                  params.xtol, support, &logger);

                    break;
                }
                case Regularizer::CONST_MAG: {
                    std::cout << "Starting constant magnitude + L2 reconstruction "
                                 "on GPU ...\n";
                    // constrain only the block that survives the crop below;
                    // the padding stays free to absorb out-of-view
                    // contributions
                    recon = opt::constmag_l2<T>(
                        A, yT, x0, params.lambda_l2, params.mu, params.maxIters,
                        params.innerIters, params.tol, params.xtol, solver_mask,
                        &logger, use_support);
                    break;
                }
                case Regularizer::COMBINED: {
                    std::cout << "Starting combined constant magnitude / TV / L2 "
                                 "reconstruction on GPU ...\n";
                    opt::RegWeights<T> w{params.const_mag, params.lambda,
                                         params.lambda_l2};
                    opt::MagOptions<T> mag;
                    mag.rho = params.rho;
                    mag.rho_growth = params.rho_growth;
                    mag.rho_max = params.rho_max;
                    mag.Ms = params.Ms;
                    mag.estimate_Ms = params.estimate_Ms;
                    mag.zero_outside = use_support;
                    recon = opt::regularized_admm<T>(
                        A, yT, x0, w, params.mu, params.maxIters, params.innerIters,
                        params.tol, params.xtol, solver_mask, &logger, mag);
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
    MBIR(const std::vector<Dataset_t<float>> &datasets, const ReconParams &params,
         const Array<float> &mask);
    template std::array<Array<double>, 3>
    MBIR(const std::vector<Dataset_t<double>> &datasets, const ReconParams &params,
         const Array<double> &mask);

} // namespace tomocam::gpu
