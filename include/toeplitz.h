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
#ifndef TOEPLITZ_H
#define TOEPLITZ_H

#include <array>
#include <complex>
#include <cstdint>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "fft.h"
#include "finufft_plan.h"
#include "padding.h"
#include "polar_grid.h"
#include "projection.h"

namespace tomocam::cpu {

    // smallest 2^a * 3^b >= n -- used to pick FFT/NUFFT-friendly padded sizes
    inline size_t next_fast_dim(size_t n) {
        size_t best = 1;
        while (best < n) best <<= 1;
        for (size_t q3 = 3; q3 < best; q3 *= 3) {
            size_t x = q3;
            while (x < n) x <<= 1;
            if (x < best) best = x;
        }
        return best;
    }

    /**
     * @brief A class for representing a single point spread function kernel
     * for vector tomography.
     * @tparam T
     */
    template <typename T>
    class PointSpreadFunction {
      private:
        using complex_t = std::complex<T>;
        dims_t dims_;
        Array<complex_t> kernel_hat_;

      public:
        [[nodiscard]] dims_t dims() const { return dims_; }
        [[nodiscard]] const Array<complex_t> &kernel() const { return kernel_hat_; }

        PointSpreadFunction() = default;

        // `plan` must already be a made (make_plan) type-1, 3D plan with
        // set_points called for `grid`; the caller owns its lifetime.
        PointSpreadFunction(const PolarGrid<T> &grid, dims_t dims,
                            const std::vector<T> &weights,
                            nufft::FinufftPlanWrapper<T> &plan)
            : dims_(dims) {

            // unit strength at every non-uniform grid point, scaled per-point
            // by the supplied weight (e.g. n_hat[i]*n_hat[j])
            auto ones = Array<complex_t>::ones(grid.dims());
            for (size_t i = 0; i < grid.dims().n1; ++i) {
                auto slc = ones.slice(i, i + 1);
                auto w = weights[i];
                std::transform(std::execution::par_unseq, slc.begin(), slc.end(),
                               slc.begin(), [w](complex_t v) { return w * v; });
            }

            // NUFFT type-1: non-uniform points -> uniform grid (backprojection
            // of the weighted unit strengths)
            Array<complex_t> nufft_out(dims_);
            int ierr = plan.execute(ones.begin(), nufft_out.begin());
            if (ierr != 0) {
                throw std::runtime_error("Error in FINUFFT execute (type-1)");
            }

            // Fourier-domain kernel: R2C FFT of the real-valued PSF
            auto psf = array::to_real(nufft_out);
            kernel_hat_ = fft::fft3_r2c(psf);
        }

        Array<T> convolve(const Array<T> &input) const {

            auto orig_dims = input.dims();

            // zero-pad (right-pad) input to match PSF size
            Array<T> padded = pad3d(input, dims_, dims_t(0, 0, 0));

            // 3D R2C forward FFT, multiply in Fourier space, 3D C2R inverse FFT
            auto output_hat = fft::fft3_r2c(padded);
            output_hat *= kernel_hat_;
            Array<T> result = fft::fft3_c2r(output_hat, dims_);

            // FFTW's R2C/C2R round trip is unnormalized by the padded volume
            T norm = static_cast<T>(dims_.n1) * static_cast<T>(dims_.n2) *
                     static_cast<T>(dims_.n3);

            // The FINUFFT type-1 output is in CMCL ordering: the PSF center
            // (k=0) sits at index dims_/2 in the padded array. The valid
            // linear-convolution result therefore starts at that same offset.
            dims_t center{dims_.n1 / 2, dims_.n2 / 2, dims_.n3 / 2};
            Array<T> cropped = crop3d(result, orig_dims, center);
            cropped /= norm;
            return cropped;
        }
    };

    /**
     * @brief Owns the 6 unique PSF kernels of the symmetric [3 x 3]
     * projection tensor for vector tomography (one per unique i<=j pair) and
     * dispatches access/convolution by (i, j).
     * @tparam T
     */
    template <typename T>
    class ToeplitzVectorOp {
      private:
        using complex_t = std::complex<T>;
        dims_t dims_;
        std::array<PointSpreadFunction<T>, 6> kernels_;

        // packs the 6 unique entries of a symmetric 3x3 tensor into [0, 6)
        static constexpr int idx(int i, int j) {
            if (i > j) std::swap(i, j);
            return i == 0 ? j : (i == 1 ? j + 2 : 5);
        }

      public:
        ToeplitzVectorOp(const PolarGrid<T> &grid, const dims_t &recon_dims) {
            dims_ = {next_fast_dim(2 * recon_dims.n1 - 1),
                     next_fast_dim(2 * recon_dims.n2 - 1),
                     next_fast_dim(2 * recon_dims.n3 - 1)};

            // one type-1 plan, reused for all 6 kernels (same grid, same
            // n_modes); scoped to this constructor so it's destroyed before
            // returning, rather than kept alive in a global cache
            {
                std::array<int64_t, 3> n_modes = {static_cast<int64_t>(dims_.n3),
                                                  static_cast<int64_t>(dims_.n2),
                                                  static_cast<int64_t>(dims_.n1)};
                nufft::FinufftPlanWrapper<T> plan;
                plan.make_plan(1, 3, n_modes, 1);
                plan.set_points(grid);

                for (size_t i = 0; i < 3; ++i) {
                    for (size_t j = i; j < 3; ++j) {
                        std::vector<T> weights(grid.dims().n1);
                        for (size_t k = 0; k < grid.dims().n1; ++k) {
                            auto n_hat = beam_dir_vector(
                                grid.angle(k), grid.gamma(k), grid.beta(k));
                            weights[k] = n_hat[i] * n_hat[j];
                        }
                        kernels_[idx(i, j)] =
                            PointSpreadFunction<T>(grid, dims_, weights, plan);
                    }
                }
            } // plan destroyed here, before the constructor returns
        }

        [[nodiscard]] const PointSpreadFunction<T> &get(int i, int j) const {
            return kernels_[idx(i, j)];
        }
    };

    // A^T A x for vector tomography via the Toeplitz/PSF trick:
    // output[i] = sum_j kernels.get(i,j).convolve(x[j]).
    template <typename T>
    std::array<Array<T>, 3> sysmat(const std::array<Array<T>, 3> &x,
                                   const ToeplitzVectorOp<T> &kernels) {
        using complex_t = std::complex<T>;
        dims_t dims = kernels.get(0, 0).dims();
        dims_t orig_dims = x[0].dims();

        // forward FFT of each (padded) input once
        std::array<Array<complex_t>, 3> xhat;
        for (size_t j = 0; j < 3; ++j) {
            Array<T> padded = pad3d(x[j], dims, dims_t(0, 0, 0));
            xhat[j] = fft::fft3_r2c(padded);
        }

        // accumulate in Fourier space: yhat[i] = sum_j kernel(i,j) * xhat[j]
        std::array<Array<complex_t>, 3> yhat;
        for (size_t i = 0; i < 3; ++i) {
            yhat[i] = xhat[0].clone();
            yhat[i] *= kernels.get(i, 0).kernel();
            for (size_t j = 1; j < 3; ++j) {
                auto term = xhat[j].clone();
                term *= kernels.get(i, j).kernel();
                yhat[i] += term;
            }
        }

        // FFTW's R2C/C2R round trip is unnormalized by the padded volume
        T fft_norm = static_cast<T>(dims.n1) * static_cast<T>(dims.n2) *
                     static_cast<T>(dims.n3);

        // match the normalization convention of the direct/NUFFT sysmat (see
        // gradient.cpp): forward/adjoint each notionally divide by the volume
        // n1*n2*n3, so A^T A divides by (n1*n2*n3)^2/(n2*n3) = n1^2*n2*n3.
        T scale = static_cast<T>(orig_dims.n1) * static_cast<T>(x[0].size());

        // see PointSpreadFunction::convolve() for why the crop starts at
        // dims/2 (FINUFFT's CMCL mode ordering).
        dims_t center{dims.n1 / 2, dims.n2 / 2, dims.n3 / 2};

        std::array<Array<T>, 3> output;
        for (size_t i = 0; i < 3; ++i) {
            Array<T> result = fft::fft3_c2r(yhat[i], dims);
            output[i] = crop3d(result, orig_dims, center);
            output[i] /= (fft_norm * scale);
        }
        return output;
    }
} // namespace tomocam::cpu

#endif // TOEPLITZ_H
