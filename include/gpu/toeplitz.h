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
#ifndef TOMOCAM_GPU_TOEPLITZ_H
#define TOMOCAM_GPU_TOEPLITZ_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <cuda/std/complex>

#include "dtypes.h"
#include "gpu/cufinufft_plan.h"
#include "gpu/device_array.h"
#include "gpu/polar_grid.h"
#include "gpu/vec_array.h"

namespace tomocam::gpu {

    template <typename T>
    using Complex = cuda::std::complex<T>;

    // smallest 2^a * 3^b >= n -- used to pick FFT/NUFFT-friendly padded sizes
    // (identical to cpu::next_fast_dim, include/toeplitz.h)
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

    // PSF-construction strategy for ToeplitzVectorOp:
    //  - SEQUENTIAL: one type-1 NUFFT execute() per kernel (ntrans=1), reusing
    //    a single scoped plan -- mirrors cpu::ToeplitzVectorOp exactly. Safer
    //    on small-GPU memory budgets; the default.
    //  - BATCHED: one type-1 NUFFT execute() with ntrans=6, producing all 6
    //    backprojected volumes at once. Fewer kernel launches, but cufinufft's
    //    batched spreading/fine-grid workspace scales with ntrans -- roughly
    //    6x the memory footprint of a single SEQUENTIAL kernel build, which
    //    risks OOM on smaller GPUs.
    enum class ToeplitzMode { SEQUENTIAL, BATCHED };

    /**
     * @brief GPU analog of cpu::PointSpreadFunction: a single PSF kernel
     * (Fourier-domain) for vector tomography, device-resident.
     * @tparam T
     */
    template <typename T>
    class PointSpreadFunction {
      private:
        dims_t dims_;
        DeviceArray<Complex<T>> kernel_hat_;

        PointSpreadFunction(dims_t dims, DeviceArray<Complex<T>> &&kernel_hat)
            : dims_(dims), kernel_hat_(std::move(kernel_hat)) {}

      public:
        PointSpreadFunction() = default;

        [[nodiscard]] dims_t dims() const { return dims_; }
        [[nodiscard]] const DeviceArray<Complex<T>> &kernel() const {
            return kernel_hat_;
        }

        // `plan` must already be a made (constructed), set_points()-called
        // type-1, 3D, ntrans=1 plan for `grid`; the caller owns its lifetime.
        // Builds unit-strength non-uniform input scaled per-row by `weights`
        // and masked by grid.w, backprojects via `plan`, and stores the
        // resulting Fourier-domain kernel.
        PointSpreadFunction(const PolarGrid<T> &grid, dims_t dims,
                            const std::vector<T> &weights,
                            nufft::cuFinfftPlanWrapper<T> &plan);

        // Internal: build kernel_hat_ directly from an already-computed
        // type-1 NUFFT backprojection (used by ToeplitzVectorOp's BATCHED
        // path, where all 6 kernels' backprojections come from one
        // ntrans=6 execute() call rather than one execute() each).
        static PointSpreadFunction<T>
        from_backprojection(dims_t dims, DeviceArray<Complex<T>> &&nufft_out);

        // sum of PSFs on the same grid dims: A^T A is linear in the kernel, so
        // datasets' kernels can be merged into one.
        PointSpreadFunction<T> &operator+=(const PointSpreadFunction<T> &rhs) {
            kernel_hat_ += rhs.kernel_hat_;
            return *this;
        }
        [[nodiscard]] bool empty() const { return kernel_hat_.size() == 0; }
    };

    /**
     * @brief GPU analog of cpu::ToeplitzVectorOp: owns the 6 unique PSF
     * kernels of the symmetric [3 x 3] projection tensor for vector
     * tomography (one per unique i<=j pair) and dispatches access by (i, j).
     * @tparam T
     */
    template <typename T>
    class ToeplitzVectorOp {
      private:
        dims_t dims_;
        std::array<PointSpreadFunction<T>, 6> kernels_;

        // packs the 6 unique entries of a symmetric 3x3 tensor into [0, 6)
        static constexpr int idx(int i, int j) {
            if (i > j) {
                int tmp = i;
                i = j;
                j = tmp;
            }
            return i == 0 ? j : (i == 1 ? j + 2 : 5);
        }

        // builds the 6 kernels of `grid` and adds them into kernels_
        void add_grid(const PolarGrid<T> &grid, ToeplitzMode mode, int device);

      public:
        // gpu_id: scaffolding for future multi-GPU dataset placement. It is
        // forwarded to the cufinufft plan(s) used to build the 6 PSF
        // kernels (matching cuFinfftPlanWrapper's existing gpu_id param),
        // but this constructor does not itself call cudaSetDevice, spawn a
        // stream, or otherwise orchestrate cross-device work -- callers
        // (e.g. gpu::mbir.cu) currently always pass -1 (use the current
        // device) for every dataset. -1 means "use the current device".
        ToeplitzVectorOp(const PolarGrid<T> &grid, const dims_t &recon_dims,
                         ToeplitzMode mode = ToeplitzMode::SEQUENTIAL,
                         int gpu_id = -1);

        // Combined operator for several grids (e.g. one per dataset): the
        // kernels are summed in place, so device memory holds 6 kernels
        // regardless of the number of grids, instead of 6 per grid.
        ToeplitzVectorOp(const std::vector<PolarGrid<T>> &grids,
                         const dims_t &recon_dims,
                         ToeplitzMode mode = ToeplitzMode::SEQUENTIAL,
                         int gpu_id = -1);

        [[nodiscard]] const PointSpreadFunction<T> &get(int i, int j) const {
            return kernels_[idx(i, j)];
        }
    };

    // A^T A x for vector tomography via the Toeplitz/PSF trick:
    // output[i] = sum_j kernels.get(i,j).kernel() (*) x[j], done in Fourier
    // space. Normalization matches cpu::sysmat (include/toeplitz.h).
    template <typename T>
    VecArray<T> sysmat(const VecArray<T> &x, const ToeplitzVectorOp<T> &kernels);

} // namespace tomocam::gpu

#endif // TOMOCAM_GPU_TOEPLITZ_H
