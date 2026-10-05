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
#include <utility>
#include <vector>

#include <cuda/std/complex>
#include <thrust/device_vector.h>

#include "projection.h" // tomocam::beam_dir_vector (host)

#include "gpu/cufinufft_plan.h"
#include "gpu/device_array.h"
#include "gpu/device_array_ops.h"
#include "gpu/device_ptr.h"
#include "gpu/fft.h"
#include "gpu/memory.h"
#include "gpu/padding.h"
#include "gpu/polar_grid.h"
#include "gpu/toeplitz.h"
#include "gpu/utils.h"
#include "gpu/vec_array.h"

namespace tomocam::gpu {

    // Scales row i (projection index) of `data` by weights[i] * mask[i,j,k],
    // mirroring cpu::PointSpreadFunction's per-row grid.rows.slice()-based
    // std::transform (include/toeplitz.h) -- GPU has no ragged-rows
    // structure, but rows are dense here (dims = nangles x nrows x ncols,
    // uniform row length), so this is a direct per-row broadcast-multiply.
    template <typename T>
    __global__ void scale_rows_kernel(DevicePtr<Complex<T>> data, const T *weights,
                                      DevicePtr<const T> mask) {
        auto idx = Index3D();
        if (idx < data.dims()) { data[idx] *= weights[idx.x] * mask[idx]; }
    }

    template <typename T>
    void scale_rows(DeviceArray<Complex<T>> &data,
                    const thrust::device_vector<T> &weights,
                    const DeviceArray<T> &mask) {
        dim3 threads(1, 16, 16);
        auto g = make_grid(data.dims(), threads);
        scale_rows_kernel<T>
            <<<g, threads>>>(data, thrust::raw_pointer_cast(weights.data()), mask);
        SAFE_CALL(cudaGetLastError());
    }

    // ---------------------------------------------------------------------
    // PointSpreadFunction
    // ---------------------------------------------------------------------

    template <typename T>
    PointSpreadFunction<T> PointSpreadFunction<T>::from_backprojection(
        dims_t dims, DeviceArray<Complex<T>> &&nufft_out) {
        auto psf = gpu::array::to_real(nufft_out);
        auto kernel_hat = gpu::fft::rfft3d(psf);
        return PointSpreadFunction<T>(dims, std::move(kernel_hat));
    }

    template <typename T>
    PointSpreadFunction<T>::PointSpreadFunction(const PolarGrid<T> &grid,
                                                dims_t dims,
                                                const std::vector<T> &weights,
                                                nufft::cuFinfftPlanWrapper<T> &plan)
        : dims_(dims) {

        // unit strength at every non-uniform grid point, scaled per-row by
        // the supplied weight (e.g. n_hat[i]*n_hat[j]) and masked by grid.w
        // so points outside [-pi,pi] don't contribute to the PSF.
        auto ones = DeviceArray<Complex<T>>::ones(grid.dims());
        thrust::device_vector<T> d_weights(weights);
        scale_rows(ones, d_weights, grid.w);

        // NUFFT type-1: non-uniform points -> uniform grid (backprojection
        // of the weighted unit strengths)
        DeviceArray<Complex<T>> nufft_out(dims_);
        plan.execute(ones.data(), nufft_out.data());

        *this = from_backprojection(dims_, std::move(nufft_out));
    }

    // ---------------------------------------------------------------------
    // ToeplitzVectorOp
    // ---------------------------------------------------------------------

    template <typename T>
    ToeplitzVectorOp<T>::ToeplitzVectorOp(const PolarGrid<T> &grid,
                                          const dims_t &recon_dims,
                                          ToeplitzMode mode, int gpu_id) {
        dims_ = {next_fast_dim(2 * recon_dims.n1 - 1),
                 next_fast_dim(2 * recon_dims.n2 - 1),
                 next_fast_dim(2 * recon_dims.n3 - 1)};

        int device = gpu_id;
        if (device < 0) { SAFE_CALL(cudaGetDevice(&device)); }
        add_grid(grid, mode, device);
    }

    template <typename T>
    ToeplitzVectorOp<T>::ToeplitzVectorOp(const std::vector<PolarGrid<T>> &grids,
                                          const dims_t &recon_dims,
                                          ToeplitzMode mode, int gpu_id) {
        dims_ = {next_fast_dim(2 * recon_dims.n1 - 1),
                 next_fast_dim(2 * recon_dims.n2 - 1),
                 next_fast_dim(2 * recon_dims.n3 - 1)};

        int device = gpu_id;
        if (device < 0) { SAFE_CALL(cudaGetDevice(&device)); }
        for (const auto &grid : grids) { add_grid(grid, mode, device); }
    }

    template <typename T>
    void ToeplitzVectorOp<T>::add_grid(const PolarGrid<T> &grid, ToeplitzMode mode,
                                       int device) {
        auto accumulate = [this](int k, PointSpreadFunction<T> &&psf) {
            if (kernels_[k].empty()) {
                kernels_[k] = std::move(psf);
            } else {
                kernels_[k] += psf;
            }
        };

        std::array<int64_t, 3> n_modes = {static_cast<int64_t>(dims_.n3),
                                          static_cast<int64_t>(dims_.n2),
                                          static_cast<int64_t>(dims_.n1)};

        // per-row weights for all 6 unique (i,j) pairs, computed on the host
        // (O(nprojs) scalar work) -- same reason cpu::ToeplitzVectorOp does
        // this host-side. Requires PolarGrid's host-resident theta/gammas/betas.
        std::array<std::vector<T>, 6> all_weights;
        for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
                std::vector<T> weights(grid.nprojs());
                for (size_t k = 0; k < grid.nprojs(); ++k) {
                    auto n_hat = tomocam::beam_dir_vector(
                        grid.angle(k), grid.gamma(k), grid.beta(k));
                    weights[k] = n_hat[i] * n_hat[j];
                }
                all_weights[idx(i, j)] = std::move(weights);
            }
        }

        if (mode == ToeplitzMode::SEQUENTIAL) {
            // one type-1 plan, reused for all 6 kernels (same grid, same
            // n_modes); scoped to this constructor so it's destroyed before
            // returning, rather than kept alive in a global cache.
            nufft::cuFinfftPlanWrapper<T> plan(1, 3, n_modes, 1, device,
                                               /*ntrans=*/1);
            plan.set_points(grid);

            for (int k = 0; k < 6; ++k) {
                accumulate(
                    k, PointSpreadFunction<T>(grid, dims_, all_weights[k], plan));
            }
        } else { // BATCHED
            // one type-1 plan with ntrans=6: builds all 6 kernels in a single
            // execute() call, sharing the non-uniform points across the
            // batch. Trades fewer launches for ~6x the batched spreading/
            // fine-grid workspace vs SEQUENTIAL -- see ToeplitzMode's doc
            // comment (include/gpu/toeplitz.h) for the memory tradeoff.
            nufft::cuFinfftPlanWrapper<T> plan(1, 3, n_modes, 1, device,
                                               /*ntrans=*/6);
            plan.set_points(grid);

            size_t npts = grid.npts;
            size_t vol_size = dims_.n1 * dims_.n2 * dims_.n3;

            auto strengths = memory::make_cunique_ptr<Complex<T>>(6 * npts);
            for (int k = 0; k < 6; ++k) {
                auto ones = DeviceArray<Complex<T>>::ones(grid.dims());
                thrust::device_vector<T> d_weights(all_weights[k]);
                scale_rows(ones, d_weights, grid.w);
                copyD2D(strengths.get() + k * npts, ones.data(),
                        npts * sizeof(Complex<T>));
            }

            auto nufft_out_all = memory::make_cunique_ptr<Complex<T>>(6 * vol_size);
            plan.execute(strengths.get(), nufft_out_all.get());

            for (int k = 0; k < 6; ++k) {
                DeviceArray<Complex<T>> block(dims_);
                copyD2D(block.data(), nufft_out_all.get() + k * vol_size,
                        vol_size * sizeof(Complex<T>));
                accumulate(k, PointSpreadFunction<T>::from_backprojection(
                                  dims_, std::move(block)));
            }
        }
    }

    // ---------------------------------------------------------------------
    // sysmat
    // ---------------------------------------------------------------------

    template <typename T>
    VecArray<T> sysmat(const VecArray<T> &x, const ToeplitzVectorOp<T> &kernels) {
        dims_t dims = kernels.get(0, 0).dims();
        dims_t orig_dims = x[0].dims();

        // forward R2C FFT of each (padded) input once
        std::array<DeviceArray<Complex<T>>, 3> xhat;
        for (size_t j = 0; j < 3; ++j) {
            auto padded = gpu::pad3d(x[j], dims, dims_t(0, 0, 0));
            xhat[j] = gpu::fft::rfft3d(padded);
        }

        // accumulate in Fourier space: yhat[i] = sum_j kernel(i,j) * xhat[j]
        std::array<DeviceArray<Complex<T>>, 3> yhat;
        for (size_t i = 0; i < 3; ++i) {
            yhat[i] = xhat[0].clone();
            yhat[i] *= kernels.get(i, 0).kernel();
            for (size_t j = 1; j < 3; ++j) {
                auto term = xhat[j].clone();
                term *= kernels.get(i, j).kernel();
                yhat[i] += term;
            }
        }

        // match cpu::sysmat's normalization convention (include/toeplitz.h):
        // forward/adjoint each notionally divide by the volume n1*n2*n3, so
        // A^T A divides by (n1*n2*n3)^2/(n2*n3) = n1^2*n2*n3, plus the padded-
        // volume factor left unnormalized by the FFT/IFFT round trip.
        T fft_norm = static_cast<T>(dims.n1) * static_cast<T>(dims.n2) *
                     static_cast<T>(dims.n3);
        T scale = static_cast<T>(orig_dims.n1) * static_cast<T>(x[0].size());

        // see PointSpreadFunction's backprojection (FINUFFT's CMCL mode
        // ordering) for why the crop starts at dims/2.
        dims_t center{dims.n1 / 2, dims.n2 / 2, dims.n3 / 2};

        VecArray<T> output;
        for (size_t i = 0; i < 3; ++i) {
            auto result = gpu::fft::irfft3d(yhat[i], dims);
            output[i] = gpu::crop3d(result, orig_dims, center);
            output[i] /= (fft_norm * scale);
        }
        return output;
    }

    // Explicit instantiations
    template class PointSpreadFunction<float>;
    template class PointSpreadFunction<double>;
    template class ToeplitzVectorOp<float>;
    template class ToeplitzVectorOp<double>;

    template VecArray<float> sysmat(const VecArray<float> &x,
                                    const ToeplitzVectorOp<float> &kernels);
    template VecArray<double> sysmat(const VecArray<double> &x,
                                     const ToeplitzVectorOp<double> &kernels);

} // namespace tomocam::gpu
