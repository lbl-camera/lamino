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

// Compare the GPU residual computation r = y - A(x) against the CPU version,
// where A(x) = sum_j [R_j^T R_j x] for 3 datasets with different tilt angles
// (gamma = 0°, 45°, -45°), matching the XMCD configuration.
//
// Both functions implement the multi-dataset normal equations operator.
// The CPU version divides by scale per dataset; the GPU version does not.
// This test will surface that difference explicitly via a reported ratio,
// in addition to computing the per-component relative L2 error.
// Pass criterion: relative L2 error < 1e-3 (after accounting for the
// scale factor, if present).

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "polar_grid.h"

// Forward-declare the CPU sysmat to avoid pulling in tomocam.h → config.h → toml++,
// which does not compile cleanly under nvcc.
namespace tomocam {
    template <typename T>
    std::array<Array<T>, 3> sysmat(const std::array<Array<T>, 3> &x,
                                   const PolarGrid<T> &grid, T gamma);
}

#include "gpu/device_array.h"
#include "gpu/polar_grid.h"
#include "gpu/projection.h"
#include "gpu/vec_array.h"

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

using namespace tomocam;

int main() {

    // -----------------------------------------------------------------------
    // Problem setup: small 3D field + 141 evenly-spaced projection angles
    // 3 datasets with gamma = 0°, 45°, -45°
    // -----------------------------------------------------------------------
    dims_t dims{21, 511, 511};
    constexpr size_t ntheta = 141;
    constexpr float DEG = static_cast<float>(M_PI) / 180.f;
    const std::array<float, 3> gammas = {0.f, 45.f * DEG, -45.f * DEG};

    std::vector<float> theta(ntheta);
    for (size_t i = 0; i < ntheta; i++)
        theta[i] = static_cast<float>(i) * 2.f * static_cast<float>(M_PI) /
                   static_cast<float>(ntheta);

    // Build CPU and GPU grids for all 3 datasets
    std::array<PolarGrid<float>, 3> cpu_grids;
    std::array<gpu::PolarGrid<float>, 3> gpu_grids;
    for (size_t j = 0; j < 3; j++) {
        cpu_grids[j] = PolarGrid<float>(theta, dims.n2, dims.n3, gammas[j]);
        gpu_grids[j] = gpu::PolarGrid<float>(theta, gammas[j], dims.n2, dims.n3);
    }

    // Random 3-component fields on CPU
    std::array<Array<float>, 3> x_cpu;
    std::array<Array<float>, 3> y_cpu;
    for (size_t i = 0; i < 3; i++) {
        x_cpu[i] = Array<float>::random(dims);
        y_cpu[i] = Array<float>::random(dims);
    }

    // Copy to GPU VecArray
    gpu::VecArray<float> x_gpu{gpu::DeviceArray<float>(x_cpu[0]),
                               gpu::DeviceArray<float>(x_cpu[1]),
                               gpu::DeviceArray<float>(x_cpu[2])};
    gpu::VecArray<float> y_gpu{gpu::DeviceArray<float>(y_cpu[0]),
                               gpu::DeviceArray<float>(y_cpu[1]),
                               gpu::DeviceArray<float>(y_cpu[2])};

    // -----------------------------------------------------------------------
    // CPU: compute r = y - A(x) where A(x) = sum_j sysmat(x, grid_j, gamma_j)
    // -----------------------------------------------------------------------
    std::cout << "Running CPU multi-dataset sysmat ...\n";
    auto t_cpu_start = std::chrono::high_resolution_clock::now();

    auto Ax_cpu = tomocam::sysmat(x_cpu, cpu_grids[0], gammas[0]);
    for (size_t j = 1; j < 3; j++) {
        auto tmp = tomocam::sysmat(x_cpu, cpu_grids[j], gammas[j]);
        for (size_t i = 0; i < 3; i++) Ax_cpu[i] += tmp[i];
    }

    std::array<Array<float>, 3> r_cpu;
    for (size_t i = 0; i < 3; i++) r_cpu[i] = y_cpu[i] - Ax_cpu[i];

    // Compute residual norm to ensure full completion
    float r_cpu_norm = 0.f;
    for (size_t i = 0; i < 3; i++) {
        float ni = array::norm2(r_cpu[i]);
        r_cpu_norm += ni * ni;
    }
    r_cpu_norm = std::sqrt(r_cpu_norm);

    auto t_cpu_end = std::chrono::high_resolution_clock::now();
    auto t_cpu_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        t_cpu_end - t_cpu_start)
                        .count();

    // -----------------------------------------------------------------------
    // GPU: compute r = y - A(x) where A(x) = sum_j sysmat(x, grid_j, gamma_j)
    // -----------------------------------------------------------------------
    std::cout << "Running GPU multi-dataset sysmat ...\n";
    auto t_gpu_start = std::chrono::high_resolution_clock::now();

    auto Ax_gpu = gpu::sysmat(x_gpu, gpu_grids[0], gammas[0]);
    for (size_t j = 1; j < 3; j++) {
        auto tmp = gpu::sysmat(x_gpu, gpu_grids[j], gammas[j]);
        Ax_gpu += tmp;
    }
    auto r_gpu = y_gpu - Ax_gpu;

    // Compute residual norm to ensure full completion and synchronization
    float r_gpu_norm = 0.f;
    for (size_t i = 0; i < 3; i++) {
        float ni = gpu::array::norm2(r_gpu[i]);
        r_gpu_norm += ni * ni;
    }
    r_gpu_norm = std::sqrt(r_gpu_norm);

#ifdef __CUDACC__
    cudaDeviceSynchronize();
#endif
    auto t_gpu_end = std::chrono::high_resolution_clock::now();
    auto t_gpu_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                        t_gpu_end - t_gpu_start)
                        .count();

    // -----------------------------------------------------------------------
    // Compare per-component residual
    // -----------------------------------------------------------------------
    std::cout << "\n--- Comparison of r = y - A(x) ---\n";
    bool all_passed = true;
    constexpr float pass_thresh = 1e-3f;

    for (size_t j = 0; j < 3; j++) {
        auto gpu_host = r_gpu[j].to_host();

        float norm_cpu = array::norm2(r_cpu[j]);
        float norm_gpu = array::norm2(gpu_host);

        // Relative difference in L2 norm (quick scale-factor check)
        float ratio = (norm_gpu > 1e-12f) ? norm_cpu / norm_gpu : 0.f;
        // steop if ration is too far from 1.0, indicating a likely bug
        if (ratio < 0.99f || ratio > 1.01f) {
            std::cout << std::format(
                "  component {}: ||cpu||={:.3e}  ||gpu||={:.3e}  ratio={:.4f}  "
                "FAILED (ratio out of bounds)\n",
                j, norm_cpu, norm_gpu, ratio);
            all_passed = false;
            continue;
        }

        // Element-wise difference after normalising GPU output to CPU scale
        // so we can tell whether the error is purely from missing normalisation
        // or from a genuine numerical mismatch.
        auto diff_raw = r_cpu[j] - gpu_host;
        float rel_err_raw = array::norm2(diff_raw) / (norm_cpu + 1e-8f);

        // Difference after rescaling GPU output by the ratio
        Array<float> gpu_scaled(dims);
        for (size_t k = 0; k < gpu_host.size(); k++)
            gpu_scaled[k] = gpu_host[k] * ratio;
        auto diff_scaled = r_cpu[j] - gpu_scaled;
        float rel_err_scaled = array::norm2(diff_scaled) / (norm_cpu + 1e-8f);

        bool passed = rel_err_scaled < pass_thresh;
        if (!passed) all_passed = false;

        std::cout << std::format(
            "  component {}: ||cpu||={:.3e}  ||gpu||={:.3e}  ratio={:.4f}"
            "  rel_err(raw)={:.3e}  rel_err(scaled)={:.3e}  {}\n",
            j, norm_cpu, norm_gpu, ratio, rel_err_raw, rel_err_scaled,
            passed ? "PASSED" : "FAILED");
    }

    // -----------------------------------------------------------------------
    // Runtime summary
    // -----------------------------------------------------------------------
    std::cout << "\n--- Timing Summary ---\n";
    std::cout << std::format("  CPU time: {:.2f} ms\n",
                             static_cast<double>(t_cpu_ms));
    std::cout << std::format("  GPU time: {:.2f} ms\n",
                             static_cast<double>(t_gpu_ms));
    if (t_gpu_ms > 0) {
        double speedup =
            static_cast<double>(t_cpu_ms) / static_cast<double>(t_gpu_ms);
        std::cout << std::format("  Speedup:  {:.2f}x\n", speedup);
    }

    if (all_passed) {
        std::cout << "\ntest_gpu_gradient: PASSED\n";
        // Diagnose normalisation discrepancy for visibility
        float norm0_cpu = array::norm2(r_cpu[0]);
        auto gpu0_host = r_gpu[0].to_host();
        float norm0_gpu = array::norm2(gpu0_host);
        if (norm0_gpu > 1e-12f) {
            float scale = static_cast<float>(cpu_grids[0].size()) /
                          static_cast<float>(cpu_grids[0].nprojs());
            std::cout << std::format("  Note: expected scale per dataset = {:.4f},  "
                                     "observed cpu/gpu ratio = {:.4f}\n",
                                     scale, norm0_cpu / norm0_gpu);
        }
    } else {
        std::cout << "\ntest_gpu_gradient: FAILED\n";
    }

    return all_passed ? 0 : 1;
}
