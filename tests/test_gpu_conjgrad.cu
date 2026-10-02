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

// Compare the GPU CG solver (src/gpu/conjgrad.cu) against the CPU CG solver
// (src/conjgrad.cpp) on an identical diagonal-matrix problem.
//
// The same diagonal operator A(x)[i] = diag[i] * x[i] is implemented on both
// CPU (via std::array<Array<T>,3>) and GPU (via thrust::transform on VecArray).
// Both solvers are initialized from the same data.  Pass criterion: max
// per-element relative error < 1e-4.

#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <vector>

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/transform.h>

#include "array.h"
#include "array_ops.h"
#include "mask.h"
#include "optimize.h"

#include "gpu/device_array.h"
#include "gpu/gpu_opt.h"
#include "gpu/vec_array.h"

using namespace tomocam;

// Build a CPU VecArray from separate CPU arrays (each cloned into the slot)
static opt::VecArray<float>
make_cpu_vec(const Array<float> &a, const Array<float> &b, const Array<float> &c) {
    opt::VecArray<float> v;
    v[0] = a.clone();
    v[1] = b.clone();
    v[2] = c.clone();
    return v;
}

int main() {

    // -----------------------------------------------------------------------
    // Problem setup: diagonal matrix A, diagMat[i] = i+1, condition number = N
    // -----------------------------------------------------------------------
    constexpr size_t N = 100;
    dims_t dims{1, 1, N};

    std::vector<float> diagMat(N);
    for (size_t i = 0; i < N; i++) diagMat[i] = static_cast<float>(i + 1);

    // Fixed true solution and RHS
    auto x_true = Array<float>::random(dims);

    Array<float> b_host(dims);
    for (size_t i = 0; i < N; i++)
        b_host[{0, 0, i}] = diagMat[i] * x_true[{0, 0, i}];

    // Initial guess (same for both CPU and GPU)
    auto x0_host = Array<float>::random(dims);

    // -----------------------------------------------------------------------
    // CPU solve
    // -----------------------------------------------------------------------
    auto A_cpu = [&](const opt::VecArray<float> &v) {
        opt::VecArray<float> Av;
        for (size_t j = 0; j < 3; j++) {
            Av[j] = Array<float>(dims);
            for (size_t i = 0; i < N; i++)
                Av[j][{0, 0, i}] = diagMat[i] * v[j][{0, 0, i}];
        }
        return Av;
    };

    auto b_cpu = make_cpu_vec(b_host, b_host, b_host);
    auto x0_cpu = make_cpu_vec(x0_host, x0_host, x0_host);

    std::cout << "--- CPU CG solver ---\n";
    // xtol=1e-30: disable step-size stopping so only residual tol drives convergence
    auto mask = mask_support<float>(dims, dims);
    auto x_sol_cpu = opt::cgsolver<float>(A_cpu, b_cpu, x0_cpu,
                                          /*max_iter=*/1000, /*tol=*/1e-8f,
                                          /*xtol=*/1e-30f, mask);

    // -----------------------------------------------------------------------
    // GPU solve (same diagMat uploaded to device)
    // -----------------------------------------------------------------------
    thrust::device_vector<float> diag_dev(diagMat);

    // Upload b and x0 to GPU
    gpu::VecArray<float> b_gpu{gpu::DeviceArray<float>(b_host),
                               gpu::DeviceArray<float>(b_host),
                               gpu::DeviceArray<float>(b_host)};
    gpu::VecArray<float> x0_gpu{gpu::DeviceArray<float>(x0_host),
                                gpu::DeviceArray<float>(x0_host),
                                gpu::DeviceArray<float>(x0_host)};

    auto A_gpu = [&](const gpu::VecArray<float> &v) {
        gpu::VecArray<float> Av;
        for (size_t j = 0; j < 3; j++) {
            gpu::DeviceArray<float> col(dims);
            thrust::transform(v[j].begin(), v[j].end(), diag_dev.begin(),
                              col.begin(), thrust::multiplies<float>());
            Av[j] = std::move(col);
        }
        return Av;
    };

    std::cout << "--- GPU CG solver ---\n";
    auto x_sol_gpu = gpu::opt::cgsolver<float>(A_gpu, b_gpu, x0_gpu,
                                               /*max_iter=*/1000, /*tol=*/1e-8f,
                                               /*xtol=*/1e-30f, dims);

    // -----------------------------------------------------------------------
    // Compare CPU vs GPU solutions
    // -----------------------------------------------------------------------
    std::cout << "\n--- Comparison ---\n";
    bool all_passed = true;
    constexpr float pass_thresh = 1e-4f;

    for (size_t j = 0; j < 3; j++) {
        auto gpu_host = x_sol_gpu[j].to_host();
        float max_rel_err = 0.f;
        for (size_t i = 0; i < N; i++) {
            float cpu_val = x_sol_cpu[j][{0, 0, i}];
            float gpu_val = gpu_host[{0, 0, i}];
            float rel = std::abs(cpu_val - gpu_val) / (std::abs(cpu_val) + 1e-8f);
            if (rel > max_rel_err) max_rel_err = rel;
        }
        bool passed = max_rel_err < pass_thresh;
        if (!passed) all_passed = false;
        std::cout << std::format("  component {}: max rel err = {:.3e}  {}\n", j,
                                 max_rel_err, passed ? "PASSED" : "FAILED");
    }

    std::cout << (all_passed ? "\ntest_gpu_conjgrad: PASSED\n"
                             : "\ntest_gpu_conjgrad: FAILED\n");
    return all_passed ? 0 : 1;
}
