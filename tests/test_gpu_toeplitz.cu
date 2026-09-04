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

// Three-way regression test for the GPU Toeplitz-trick gradient (Phase B):
//
// 1. GPU Toeplitz sysmat vs GPU direct sysmat -- both now share the same
//    normalization convention (src/gpu/gradient.cu was aligned to match
//    include/toeplitz.h's convention as part of this work), so they should
//    agree numerically, not just in direction.
// 2. GPU Toeplitz sysmat vs CPU Toeplitz sysmat -- same convention, primary
//    regression gate for the GPU port itself.
// 3. GPU Toeplitz SEQUENTIAL vs BATCHED PSF-construction modes -- same grid,
//    should produce numerically identical sysmat output; catches bugs in the
//    batched strength/output striding logic.

#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <random>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "dtypes.h"
#include "polar_grid.h"
#include "projection.h"
#include "toeplitz.h"

#include "gpu/cufft_plan_cache.h"
#include "gpu/cufinufft_plan_cache.h"
#include "gpu/device_array.h"
#include "gpu/polar_grid.h"
#include "gpu/projection.h"
#include "gpu/toeplitz.h"
#include "gpu/vec_array.h"

using namespace tomocam;

namespace {

    using Vec3 = std::array<Array<float>, 3>;

    Vec3 random_vec3(const dims_t &dims, std::mt19937 &gen) {
        std::uniform_real_distribution<float> dis(0.f, 1.f);
        Vec3 v;
        for (size_t i = 0; i < 3; ++i) {
            v[i] = Array<float>(dims);
            for (size_t k = 0; k < v[i].size(); ++k) v[i][k] = dis(gen);
        }
        return v;
    }

    Vec3 to_host(const gpu::VecArray<float> &v) {
        return {v[0].to_host(), v[1].to_host(), v[2].to_host()};
    }

    float vec_dot(const Vec3 &a, const Vec3 &b) {
        float s = 0;
        for (size_t i = 0; i < 3; ++i) s += array::dot<float>(a[i], b[i]);
        return s;
    }

    float vec_norm(const Vec3 &a) { return std::sqrt(vec_dot(a, a)); }

    float cosine_sim(const Vec3 &a, const Vec3 &b) {
        float na = vec_norm(a), nb = vec_norm(b);
        return (na > 0 && nb > 0) ? vec_dot(a, b) / (na * nb) : 0.f;
    }

    // ||a|| / ||b||
    float norm_ratio(const Vec3 &a, const Vec3 &b) {
        float nb = vec_norm(b);
        return (nb > 0) ? vec_norm(a) / nb : 0.f;
    }

    bool check(const char *label, const Vec3 &a, const Vec3 &b, float cos_thresh,
              float ratio_lo, float ratio_hi) {
        float cos = cosine_sim(a, b);
        float ratio = norm_ratio(a, b);
        bool passed = cos > cos_thresh && ratio > ratio_lo && ratio < ratio_hi;
        std::cout << std::format(
            "{}: cosine_sim={:.6f} norm_ratio(a/b)={:.4f}  {}\n", label, cos, ratio,
            passed ? "PASSED" : "FAILED");
        return passed;
    }

} // namespace

int main() {
    bool all_passed = true;

    // Small synthetic problem: odd angle/detector counts (required for
    // Hermitian-symmetric NUFFT output, see make_polar_grid_kernel).
    dims_t dims{9, 9, 9};
    constexpr size_t ntheta = 17;
    constexpr float DEG = static_cast<float>(M_PI) / 180.f;
    float gamma = 20.f * DEG;
    float beta = 15.f * DEG;

    std::vector<float> theta(ntheta);
    for (size_t i = 0; i < ntheta; ++i)
        theta[i] = static_cast<float>(i) * 2.f * static_cast<float>(M_PI) /
                   static_cast<float>(ntheta);

    PolarGrid<float> cpu_grid(theta, dims.n2, dims.n3, gamma, beta);
    gpu::PolarGrid<float> gpu_grid(theta, dims.n2, dims.n3, gamma, beta);

    std::cout << "Building CPU Toeplitz kernels ...\n";
    cpu::ToeplitzVectorOp<float> cpu_kernels(cpu_grid, dims);

    std::cout << "Building GPU Toeplitz kernels (SEQUENTIAL) ...\n";
    gpu::ToeplitzVectorOp<float> gpu_kernels_seq(gpu_grid, dims,
                                                 gpu::ToeplitzMode::SEQUENTIAL);

    std::cout << "Building GPU Toeplitz kernels (BATCHED) ...\n";
    gpu::ToeplitzVectorOp<float> gpu_kernels_batched(gpu_grid, dims,
                                                     gpu::ToeplitzMode::BATCHED);

    std::mt19937 gen(42);
    Vec3 x_cpu = random_vec3(dims, gen);
    gpu::VecArray<float> x_gpu{gpu::DeviceArray<float>(x_cpu[0]),
                               gpu::DeviceArray<float>(x_cpu[1]),
                               gpu::DeviceArray<float>(x_cpu[2])};

    auto Ax_cpu_toeplitz = cpu::sysmat(x_cpu, cpu_kernels);
    auto Ax_cpu_direct = tomocam::sysmat(x_cpu, cpu_grid);
    auto Ax_gpu_toeplitz_seq = to_host(gpu::sysmat(x_gpu, gpu_kernels_seq));
    auto Ax_gpu_toeplitz_batched = to_host(gpu::sysmat(x_gpu, gpu_kernels_batched));
    auto Ax_gpu_direct = to_host(gpu::sysmat(x_gpu, gpu_grid));

    std::cout << "\n--- Sanity: CPU Toeplitz vs CPU direct (reference) ---\n";
    check("CPU Toeplitz vs CPU direct", Ax_cpu_toeplitz, Ax_cpu_direct, 0.999f,
         0.99f, 1.01f);

    std::cout << "\n--- Test 1: GPU Toeplitz vs GPU direct ---\n";
    if (!check("GPU Toeplitz(seq) vs GPU direct", Ax_gpu_toeplitz_seq, Ax_gpu_direct,
              0.999f, 0.99f, 1.01f))
        all_passed = false;

    std::cout << "\n--- Test 2: GPU Toeplitz vs CPU Toeplitz "
                 "(primary regression gate) ---\n";
    if (!check("GPU Toeplitz(seq) vs CPU Toeplitz", Ax_gpu_toeplitz_seq,
              Ax_cpu_toeplitz, 0.999f, 0.99f, 1.01f))
        all_passed = false;

    std::cout << "\n--- Test 3: SEQUENTIAL vs BATCHED PSF construction ---\n";
    if (!check("GPU Toeplitz(seq) vs GPU Toeplitz(batched)", Ax_gpu_toeplitz_seq,
              Ax_gpu_toeplitz_batched, 0.999f, 0.99f, 1.01f))
        all_passed = false;

    std::cout << (all_passed ? "\ntest_gpu_toeplitz: PASSED\n"
                             : "\ntest_gpu_toeplitz: FAILED\n");
    std::cout.flush();

    // Tear down plan caches explicitly while the CUDA context is still fully
    // alive -- their global-static destructors otherwise run during process
    // exit / driver shutdown and abort (see tomocam::gpu::MBIR for the same
    // pattern in src/gpu/mbir.cu).
    tomocam::gpu::nufft::plans::cache<float>.clear();
    tomocam::gpu::fft::plans::cache<float>.clear();

    return all_passed ? 0 : 1;
}
