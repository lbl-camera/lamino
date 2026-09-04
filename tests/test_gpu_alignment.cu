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

// Regression + parity test for the GPU alignment-parameter work (Phase A):
// gpu::PolarGrid now carries per-projection gamma/beta and builds its
// non-uniform coordinates via the same RotationTranspose math as CPU, and
// gpu::adjoint supports the same per-projection center-of-rotation phase
// shift as CPU. This test compares GPU output against CPU (the ground
// truth), not against pre-change GPU output, since the rotation kernel
// itself changed.
//
// Part 1: beta = 0, no shifts -- polar grid coordinates (x,y,z,w) must match
//         CPU exactly (regression gate for the rewritten rotation kernel).
// Part 2: nonzero gamma, beta, and populated shifts -- forward/adjoint/sysmat
//         output must match CPU within numerical tolerance.

#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "polar_grid.h"

// Forward-declare the CPU functions under test to avoid pulling in
// tomocam.h -> config.h -> toml++, which does not compile cleanly under nvcc
// (same workaround as test_gpu_gradient.cu).
namespace tomocam {
    template <typename T>
    Array<T> forward(const std::array<Array<T>, 3> &magnetization,
                     const PolarGrid<T> &pg, T gamma, T beta);

    template <typename T>
    std::array<Array<T>, 3>
    adjoint(const Array<T> &proj, const PolarGrid<T> &pg, const dims_t &recon_dims,
           T gamma, T beta, const std::vector<std::array<T, 2>> &shifts);

    template <typename T>
    std::array<Array<T>, 3> sysmat(const std::array<Array<T>, 3> &x,
                                   const PolarGrid<T> &grid, T gamma, T beta);
}

#include "gpu/cufft_plan_cache.h"
#include "gpu/cufinufft_plan_cache.h"
#include "gpu/device_array.h"
#include "gpu/polar_grid.h"
#include "gpu/projection.h"
#include "gpu/vec_array.h"

#include <thrust/device_vector.h>

using namespace tomocam;

namespace {

    template <typename T>
    T max_abs_diff(const Array<T> &a, const Array<T> &b) {
        T worst = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            worst = std::max(worst, std::abs(a[i] - b[i]));
        }
        return worst;
    }

    // GPU forward/adjoint/sysmat now match CPU's normalization convention
    // (src/gpu/projection.cu, src/gpu/gradient.cu), so both direction
    // (cosine similarity) and magnitude (norm ratio) should agree.
    template <typename T>
    T cosine_similarity(const Array<T> &cpu, const Array<T> &gpu) {
        T dot = 0, ncpu = 0, ngpu = 0;
        for (size_t i = 0; i < cpu.size(); ++i) {
            dot += cpu[i] * gpu[i];
            ncpu += cpu[i] * cpu[i];
            ngpu += gpu[i] * gpu[i];
        }
        return dot / (std::sqrt(ncpu) * std::sqrt(ngpu) + T(1e-12));
    }

    template <typename T>
    T norm_ratio(const Array<T> &cpu, const Array<T> &gpu) {
        T ncpu = 0, ngpu = 0;
        for (size_t i = 0; i < cpu.size(); ++i) {
            ncpu += cpu[i] * cpu[i];
            ngpu += gpu[i] * gpu[i];
        }
        return std::sqrt(ncpu) / (std::sqrt(ngpu) + T(1e-12));
    }

} // namespace

int main() {
    using T = float;
    constexpr T DEG = static_cast<T>(M_PI) / T(180);
    constexpr T coord_thresh = T(1e-4);
    constexpr T cosine_thresh = T(0.999);
    constexpr T ratio_lo = T(0.99);
    constexpr T ratio_hi = T(1.01);
    bool all_passed = true;

    // Small synthetic problem: odd angle/detector counts (required for
    // Hermitian-symmetric NUFFT output, see make_polar_grid_kernel).
    dims_t dims{9, 9, 9};
    constexpr size_t ntheta = 17;
    std::vector<T> theta(ntheta);
    for (size_t i = 0; i < ntheta; ++i)
        theta[i] = static_cast<T>(i) * T(2) * static_cast<T>(M_PI) /
                  static_cast<T>(ntheta);

    // -------------------------------------------------------------------
    // Part 1: beta = 0, gamma = 0 -- polar grid coordinates must match CPU
    // exactly (this is the regression gate for the rewritten GPU rotation
    // kernel, which now differs from pre-change GPU output).
    // -------------------------------------------------------------------
    {
        PolarGrid<T> cpu_grid(theta, dims.n2, dims.n3, T(0), T(0));
        gpu::PolarGrid<T> gpu_grid(theta, dims.n2, dims.n3, T(0), T(0));

        auto gx = gpu_grid.x.to_host();
        auto gy = gpu_grid.y.to_host();
        auto gz = gpu_grid.z.to_host();
        auto gw = gpu_grid.w.to_host();

        T dx = max_abs_diff(cpu_grid.x, gx);
        T dy = max_abs_diff(cpu_grid.y, gy);
        T dz = max_abs_diff(cpu_grid.z, gz);
        T dw = max_abs_diff(cpu_grid.w, gw);

        bool passed =
            dx < coord_thresh && dy < coord_thresh && dz < coord_thresh && dw == T(0);
        if (!passed) all_passed = false;
        std::cout << std::format(
            "Part 1 (beta=0 grid coords): max|dx|={:.3e} max|dy|={:.3e} "
            "max|dz|={:.3e} max|dw|={:.3e}  {}\n",
            dx, dy, dz, dw, passed ? "PASSED" : "FAILED");
    }

    // -------------------------------------------------------------------
    // Part 2: nonzero gamma, beta, and populated shifts -- forward/adjoint/
    // sysmat must match CPU within numerical tolerance.
    // -------------------------------------------------------------------
    {
        T gamma = T(20) * DEG;
        T beta = T(15) * DEG;

        PolarGrid<T> cpu_grid(theta, dims.n2, dims.n3, gamma, beta);
        gpu::PolarGrid<T> gpu_grid(theta, dims.n2, dims.n3, gamma, beta);

        // grid coordinates should still agree at nonzero beta
        {
            auto gx = gpu_grid.x.to_host();
            auto gy = gpu_grid.y.to_host();
            auto gz = gpu_grid.z.to_host();
            T dx = max_abs_diff(cpu_grid.x, gx);
            T dy = max_abs_diff(cpu_grid.y, gy);
            T dz = max_abs_diff(cpu_grid.z, gz);
            bool passed = dx < coord_thresh && dy < coord_thresh && dz < coord_thresh;
            if (!passed) all_passed = false;
            std::cout << std::format(
                "Part 2 (beta!=0 grid coords): max|dx|={:.3e} max|dy|={:.3e} "
                "max|dz|={:.3e}  {}\n",
                dx, dy, dz, passed ? "PASSED" : "FAILED");
        }

        // per-projection COR shifts, a couple of nonzero entries
        std::vector<std::array<T, 2>> shifts(ntheta, {T(0), T(0)});
        shifts[0] = {T(0.7), T(-0.3)};
        shifts[5] = {T(-1.1), T(0.4)};
        std::vector<T> hx(ntheta), hy(ntheta);
        for (size_t i = 0; i < ntheta; ++i) {
            hx[i] = shifts[i][0];
            hy[i] = shifts[i][1];
        }
        thrust::device_vector<T> shift_dx(hx), shift_dy(hy);

        // random 3-component field / projections
        std::array<Array<T>, 3> m_cpu;
        for (size_t i = 0; i < 3; ++i) m_cpu[i] = Array<T>::random(dims);
        gpu::VecArray<T> m_gpu{gpu::DeviceArray<T>(m_cpu[0]),
                               gpu::DeviceArray<T>(m_cpu[1]),
                               gpu::DeviceArray<T>(m_cpu[2])};

        // forward
        {
            auto proj_cpu = tomocam::forward(m_cpu, cpu_grid, gamma, beta);
            auto proj_gpu = gpu::forward(m_gpu, gpu_grid).to_host();
            T cos = cosine_similarity(proj_cpu, proj_gpu);
            T ratio = norm_ratio(proj_cpu, proj_gpu);
            bool passed = cos > cosine_thresh && ratio > ratio_lo && ratio < ratio_hi;
            if (!passed) all_passed = false;
            std::cout << std::format(
                "Part 2 (forward): cosine_sim={:.6f} norm_ratio(cpu/gpu)={:.4f}  {}\n",
                cos, ratio, passed ? "PASSED" : "FAILED");
        }

        // adjoint (exercises the phase-shift path)
        auto proj_cpu = tomocam::forward(m_cpu, cpu_grid, gamma, beta);
        auto proj_gpu_dev = gpu::forward(m_gpu, gpu_grid);
        {
            auto adj_cpu =
                tomocam::adjoint(proj_cpu, cpu_grid, dims, gamma, beta, shifts);
            auto adj_gpu = gpu::adjoint(proj_gpu_dev, gpu_grid, dims, shift_dx,
                                        shift_dy);
            bool passed = true;
            for (size_t i = 0; i < 3; ++i) {
                auto gpu_host = adj_gpu[i].to_host();
                T cos = cosine_similarity(adj_cpu[i], gpu_host);
                T ratio = norm_ratio(adj_cpu[i], gpu_host);
                bool ok =
                    cos > cosine_thresh && ratio > ratio_lo && ratio < ratio_hi;
                if (!ok) passed = false;
                std::cout << std::format(
                    "Part 2 (adjoint[{}]): cosine_sim={:.6f} "
                    "norm_ratio(cpu/gpu)={:.4f}  {}\n",
                    i, cos, ratio, ok ? "PASSED" : "FAILED");
            }
            if (!passed) all_passed = false;
        }

        // sysmat
        {
            auto Ax_cpu = tomocam::sysmat(m_cpu, cpu_grid, gamma, beta);
            auto Ax_gpu = gpu::sysmat(m_gpu, gpu_grid);
            bool passed = true;
            for (size_t i = 0; i < 3; ++i) {
                auto gpu_host = Ax_gpu[i].to_host();
                T cos = cosine_similarity(Ax_cpu[i], gpu_host);
                T ratio = norm_ratio(Ax_cpu[i], gpu_host);
                bool ok =
                    cos > cosine_thresh && ratio > ratio_lo && ratio < ratio_hi;
                if (!ok) passed = false;
                std::cout << std::format(
                    "Part 2 (sysmat[{}]): cosine_sim={:.6f} "
                    "norm_ratio(cpu/gpu)={:.4f}  {}\n",
                    i, cos, ratio, ok ? "PASSED" : "FAILED");
            }
            if (!passed) all_passed = false;
        }
    }

    std::cout << (all_passed ? "\ntest_gpu_alignment: PASSED\n"
                             : "\ntest_gpu_alignment: FAILED\n");
    std::cout.flush();

    // Tear down plan caches explicitly while the CUDA context is still fully
    // alive -- their global-static destructors otherwise run during process
    // exit / driver shutdown and abort (see tomocam::gpu::MBIR for the same
    // pattern in src/gpu/mbir.cu).
    tomocam::gpu::nufft::plans::cache<float>.clear();
    tomocam::gpu::fft::plans::cache<float>.clear();

    return all_passed ? 0 : 1;
}
