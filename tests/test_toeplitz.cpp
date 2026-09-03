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

// Compare the two ways of computing A^T A x for vector (XMCD) tomography:
//   sysmat(x, PolarGrid)         -- direct NUFFT path (projection.h / gradient.cpp)
//   sysmat(x, ToeplitzVectorOp)      -- Toeplitz (PSF convolution) path (toeplitz.h)
//
// Both compute the same normal operator for the same geometry; they should
// agree in direction (cosine similarity) and, after accounting for the
// differing absolute normalization of the two paths, in shape.

#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "array.h"
#include "array_ops.h"
#include "dtypes.h"
#include "projection.h"
#include "timer.h"
#include "toeplitz.h"

using namespace tomocam;

// Array<T>::random() shares one RNG across std::execution::par_unseq
// worker threads (unrelated pre-existing data race) -- fill sequentially here
// instead.
static Array<double> random_array(const dims_t &dims, std::mt19937 &gen) {
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    Array<double> a(dims);
    for (size_t i = 0; i < a.size(); ++i) a[i] = dis(gen);
    return a;
}

static double vec_dot(const std::array<Array<double>, 3> &a,
                      const std::array<Array<double>, 3> &b) {
    double s = 0;
    for (size_t i = 0; i < 3; ++i) s += array::dot<double>(a[i], b[i]);
    return s;
}

static double vec_norm(const std::array<Array<double>, 3> &a) {
    return std::sqrt(vec_dot(a, a));
}

// Returns cosine similarity: <a,b> / (||a|| ||b||).  1.0 = identical direction.
static double cosine_sim(const std::array<Array<double>, 3> &a,
                         const std::array<Array<double>, 3> &b) {
    double na = vec_norm(a), nb = vec_norm(b);
    return (na > 0 && nb > 0) ? vec_dot(a, b) / (na * nb) : 0.0;
}

static bool test_sysmat_agreement(const dims_t &vol_dims,
                                  const std::vector<double> &theta, double gamma,
                                  double beta, double tol) {
    std::cout << std::format(
        "vol({},{},{}) nangles={:<3}  sysmat NUFFT vs Toeplitz\n", vol_dims.n1,
        vol_dims.n2, vol_dims.n3, theta.size());

    auto pg = PolarGrid<double>(theta, vol_dims.n2, vol_dims.n3, gamma, beta);

    Timer t_build;
    t_build.start();
    cpu::ToeplitzVectorOp<double> kernels(pg, vol_dims);
    t_build.stop();

    std::mt19937 gen(42);
    std::array<Array<double>, 3> x;
    for (size_t i = 0; i < 3; ++i) x[i] = random_array(vol_dims, gen);

    Timer t_nufft;
    t_nufft.start();
    auto nufft_result = sysmat(x, pg, gamma, beta);
    t_nufft.stop();

    Timer t_toeplitz;
    t_toeplitz.start();
    auto toeplitz_result = cpu::sysmat(x, kernels);
    t_toeplitz.stop();

    double speedup = (t_toeplitz.seconds() > 0)
                         ? t_nufft.seconds() / t_toeplitz.seconds()
                         : 0.0;
    std::cout << std::format(
        "  time: kernel_build={:.3f}s  nufft_sysmat={:.3f}s  "
        "toeplitz_sysmat={:.3f}s  (toeplitz apply {:.2f}x direct)\n",
        t_build.seconds(), t_nufft.seconds(), t_toeplitz.seconds(), speedup);

    double n_nufft = vec_norm(nufft_result);
    double n_toeplitz = vec_norm(toeplitz_result);

    double scale_ratio =
        vec_dot(nufft_result, toeplitz_result) / (n_toeplitz * n_toeplitz);
    double cos_sim = cosine_sim(nufft_result, toeplitz_result);

    // normalised relative difference: compare nufft with rescaled toeplitz
    std::array<Array<double>, 3> diff;
    for (size_t i = 0; i < 3; ++i) {
        auto toeplitz_scaled = toeplitz_result[i] * scale_ratio;
        diff[i] = nufft_result[i] - toeplitz_scaled;
    }
    double rel_diff_norm = vec_norm(diff) / n_nufft;

    bool shape_ok = rel_diff_norm < tol;
    bool scale_ok = std::abs(scale_ratio - 1.0) < tol;
    bool ok = shape_ok && scale_ok;

    std::cout << std::format("  ||nufft||={:.4e}  ||toeplitz||={:.4e}  "
                             "scale_ratio={:.4e}  cos_sim={:.6f}\n",
                             n_nufft, n_toeplitz, scale_ratio, cos_sim);
    std::cout << std::format("  normalised rel_diff={:.2e}  shape={} scale={}\n",
                             rel_diff_norm, shape_ok ? "OK" : "FAIL",
                             scale_ok ? "OK" : "FAIL");
    return ok;
}

int main(int argc, char **argv) {
    std::cout << "\n====== sysmat: NUFFT vs Toeplitz agreement tests (XMCD) "
                 "======\n\n";

    // n2 and n3 must be ODD -- see test_adjoint.cpp for why: the polar grid's
    // half-pixel frequency offset only gives a Hermitian-symmetric NUFFT
    // output (and thus a clean real-valued PSF) for odd n.
    const double gamma = 0.0;
    const double beta = 0.0;
    const double tol = 1e-3;
    int passed = 0, failed = 0;
    auto record = [&](bool ok) { ok ? ++passed : ++failed; };

    // Small volume, few angles
    {
        const dims_t vol_dims{5, 15, 15};
        const size_t nangles = 7;
        std::vector<double> theta(nangles);
        for (size_t i = 0; i < nangles; ++i)
            theta[i] = (static_cast<double>(i) - nangles / 2.0) * M_PI / nangles;
        record(test_sysmat_agreement(vol_dims, theta, gamma, beta, tol));
    }

    // Larger volume, more angles
    {
        const dims_t vol_dims{8, 31, 31};
        const size_t nangles = 13;
        std::vector<double> theta(nangles);
        for (size_t i = 0; i < nangles; ++i)
            theta[i] = (static_cast<double>(i) - nangles / 2.0) * M_PI / nangles;
        record(test_sysmat_agreement(vol_dims, theta, gamma, beta, tol));
    }

    // Non-zero gamma
    {
        const dims_t vol_dims{5, 15, 15};
        const size_t nangles = 9;
        const double gam = 0.3;
        std::vector<double> theta(nangles);
        for (size_t i = 0; i < nangles; ++i)
            theta[i] = (static_cast<double>(i) - nangles / 2.0) * M_PI / nangles;
        record(test_sysmat_agreement(vol_dims, theta, gam, beta, tol));
    }

    std::cout << "\nResults: " << passed << " passed, " << failed << " failed\n";

    // Optional problem-size timing comparison, e.g.:
    //   ./test_toeplitz 21 512 512
    // n2/n3 need not be odd here -- this block is timing only, not
    // tolerance-gated (see the odd-n comment above).
    if (argc > 1) {
        if (argc != 4) {
            std::cerr << "usage: " << argv[0] << " [n1 n2 n3]\n";
            return 1;
        }
        const dims_t vol_dims{static_cast<size_t>(std::stoul(argv[1])),
                              static_cast<size_t>(std::stoul(argv[2])),
                              static_cast<size_t>(std::stoul(argv[3]))};
        const size_t nangles = 91;
        std::cout << std::format("\n--- Large problem size timing ({},{},{}) ---\n",
                                 vol_dims.n1, vol_dims.n2, vol_dims.n3);
        std::vector<double> theta(nangles);
        for (size_t i = 0; i < nangles; ++i)
            theta[i] = (static_cast<double>(i) - nangles / 2.0) * M_PI / nangles;
        test_sysmat_agreement(vol_dims, theta, gamma, beta, tol);
    }

    return (failed == 0) ? 0 : 1;
}
